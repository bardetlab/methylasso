#include <Eigen/Core>
#include <RcppEigen.h>

#include "offset.hpp"
#include "binner_helpers.hpp" //make_datasets_binner
#include "Observations.hpp"
#include "Design.hpp"
#include "design_helpers.hpp" //make_group_design

namespace MethyLasso {
namespace estimator {

template<>
std::shared_ptr<binned::Binner> make_binner_ptr<Offset,Mean>(const data::Observations& data, const Config<Offset,Mean>& conf) {
  //distance and bounds
  Eigen::ArrayXd pos = data.get_position();
  Eigen::ArrayXi dset = data.get_dataset();
  //place each dataset in one bin
  auto design = binned::make_even_bins(pos.minCoeff(),pos.maxCoeff(), 1);
  //build binner
  auto pos_by_dset = binned::split_by_dataset(dset, pos);
  std::shared_ptr<binned::Binner> binner_ptr = 
    binned::make_datasets_binner(pos_by_dset, &binned::make_binner_using_bin_design, design, true );  //false to keep all bins
  return(binner_ptr);
}

template<>
params::Design make_design<Offset,Mean>(const binned::Binner& binner, const Config<Offset,Mean>&,
                                        const data::DataDesign&) {
  auto design = params::make_identity_design(binner);
  return design;
}

} //namespace estimator
}

