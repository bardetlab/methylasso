#include <Eigen/Sparse>
#include <Rcpp.h>
#include <unordered_map>

#include "data_helpers.hpp"

#include "Observations.hpp"
#include "DataDesign.hpp"

namespace MethyLasso {
namespace data {

Eigen::ArrayXi group_datasets(const DataDesign& data_design, double tol) {
  //check input data
  Rcpp::Rcout << "group_datasets\n";
  if (data_design.get_num_estimators() != 1) throw std::invalid_argument("data_design must have only one column!");
  Eigen::SparseVector<double> col = data_design.get_col(0);
  //form groups
  std::unordered_map<int, std::vector<unsigned> > groups;
  for (Eigen::SparseVector<double>::InnerIterator it(col); it; ++it) {
    int coef = it.value() / tol; //cast to int so the map access does not rely on float keys
    unsigned dset_id = it.index() + 1;
    groups[coef].push_back(dset_id);
    Rcpp::Rcout << "coef=" << coef << " now adds dset=" << dset_id << "\n";
  }
  //write group number to array
  Eigen::ArrayXi group_ids(Eigen::ArrayXi::Zero(col.rows()));
  unsigned group_id = 0;
  for (auto it = groups.cbegin(); it != groups.cend(); ++it) {
    group_id++;
    for (auto d : it->second) group_ids(d-1) = group_id;
  }
  return group_ids;
}

DataDesign get_group_design(const DataDesign& data_design, const Eigen::ArrayXi& group_ids) {
  Eigen::VectorXd col(data_design.get_col(0));
  unsigned ngroups = group_ids.maxCoeff();
  Eigen::VectorXd design(Eigen::VectorXd::Zero(ngroups));
  for (unsigned i=0; i<group_ids.size(); ++i)
    if (group_ids(i)>0) //groups with ID 0 are dropped
      design(group_ids(i)-1) = col(i); //groups are formed based on equal values for col(i), so we can replace safely
  return DataDesign{design.sparseView()};    
}

}
}
