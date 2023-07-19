#include <Eigen/Core>
#include <RcppEigen.h>

#include "unique_positions.hpp"
#include "binner_helpers.hpp" //make_datasets_binner
#include "Observations.hpp"
#include "Design.hpp"
#include "design_helpers.hpp" //make_group_design
#include "DataDesign.hpp"
#include "data_helpers.hpp"

namespace MethyLasso {
namespace estimator {

template<>
std::shared_ptr<binned::Binner> make_binner_ptr<UniquePositions,Lasso>(const data::Observations& data,
                                                                       const Config<UniquePositions,Lasso>& conf) {
  //bin each dataset separately (block diagonal binner matrix)
  Eigen::ArrayXi dset = data.get_dataset();
  //binner matrix
  Eigen::ArrayXd pos = data.get_position();
  double genome_sz = (pos.maxCoeff()-pos.minCoeff()+1);
  unsigned nbins = genome_sz;
  //compute evenly spaced bins for all data
  auto design = binned::make_unique_bins(pos.cast<int>());
  //build binner
  auto pos_by_dset = binned::split_by_dataset(dset, pos);
  std::shared_ptr<binned::Binner> binner_ptr = 
    binned::make_datasets_binner(pos_by_dset, &binned::make_binner_using_bin_design, design, true );  //false to keep all bins
  return(binner_ptr);
}

template<>
params::Design make_design<UniquePositions,Lasso>(const binned::Binner& binner, const Config<UniquePositions,Lasso>& conf,
                                               const data::DataDesign& data_design) {
  auto design = params::make_group_design(binner, data_design);
  return design;
}

} //namespace estimator
}

