#include <Eigen/Sparse>
#include <Rcpp.h>
#include <map>

#include "data_design_helpers.hpp"

namespace MethyLasso {
namespace data {

DataDesign build_signal_design(const Observations& obs) {
  const auto dset = obs.get_dataset();
  const auto cond = obs.get_condition();
  //build triplet list
  std::vector<Eigen::Triplet<double> > tripletList;
  int prev_dset = 0; //datasets start at 1
  for (unsigned i=0; i<obs.get_Ndata(); ++i) {
    if (prev_dset != dset(i)) {
      prev_dset = dset(i);
      tripletList.emplace_back(dset(i)-1,cond(i)-1,1.);
    }
  }
  //create sparse matrix and return
  const unsigned num_dsets = dset.tail<1>()(0);
  const unsigned num_conditions = cond.tail<1>()(0);
  Eigen::SparseMatrix<double> ret(num_dsets,num_conditions);
  ret.setFromTriplets(tripletList.begin(), tripletList.end());
  return DataDesign(ret);
}

DataDesign build_difference_design(const Observations& obs, unsigned ref) {
  //start with signal design
  Eigen::MatrixXd design = build_signal_design(obs).get_X();
  const unsigned num_dsets = design.rows();
  const unsigned num_conditions = design.cols();
  if (num_dsets < ref || ref == 0) throw std::invalid_argument("Invalid reference dataset!");
  //remove column with reference condition and add vector of ones in the startning
  if (ref>1) design.block(0,1,num_dsets,ref-1) = design.leftCols(ref-1).eval();
  design.col(0) = Eigen::VectorXd::Ones(num_dsets);
  return DataDesign(design.sparseView());
}

}
}
