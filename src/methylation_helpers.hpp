#ifndef METHYLATION_HELPERS_HPP
#define METHYLATION_HELPERS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

#include "Estimator.hpp"
#include "EstimatorGroup.hpp"

namespace MethyLasso {

template<typename Leg, typename Method>
Rcpp::DataFrame get_estimate_dataframe(const estimator::Estimator<Leg,Method>& est) {
  const Eigen::ArrayXd support = est.get_design().get_support();
  const Eigen::ArrayXd beta = est.get_params().get_beta();
  const Eigen::ArrayXd betahat = est.get_pseudodata().get_betahat();
  const Eigen::ArrayXd pseudoweight = est.get_pseudodata().get_pseudoweight();
  const Eigen::ArrayXi one = Eigen::ArrayXi::Ones(support.size());
  return Rcpp::DataFrame::create(_["estimate_id"]=one, _["start"]=support.cast<int>(), _["beta"]=beta,
                                 _["betahat"]=betahat, _["pseudoweight"]=pseudoweight);
}

template<typename Leg, typename Method>
Rcpp::DataFrame get_estimate_dataframe(const estimator::EstimatorGroup<Leg,Method>& eg) {
  //call get_estimate_dataframe on each estimator
  std::vector<Rcpp::DataFrame> dfs;
  unsigned total_sz = 0;
  for (auto&& est : eg.get_estimators()) {
    dfs.push_back(get_estimate_dataframe(est));
    total_sz += dfs.back().nrows();
  }
  //merge in a single dataframe
  Eigen::ArrayXd support(total_sz), beta(total_sz), betahat(total_sz), pseudoweight(total_sz);
  Eigen::ArrayXi idx(total_sz);
  unsigned this_idx = 0;
  for (unsigned i=0; i<dfs.size(); ++i) {
    auto df = dfs[i];
    double this_sz = df.nrows();
    idx.segment(this_idx, this_sz) = Eigen::ArrayXi::Constant(this_sz,i+1);
    support.segment(this_idx, this_sz) = Rcpp::as<Eigen::ArrayXd>(df["start"]);
    beta.segment(this_idx, this_sz) = Rcpp::as<Eigen::ArrayXd>(df["beta"]);
    betahat.segment(this_idx, this_sz) = Rcpp::as<Eigen::ArrayXd>(df["betahat"]);
    pseudoweight.segment(this_idx, this_sz) = Rcpp::as<Eigen::ArrayXd>(df["pseudoweight"]);
    this_idx += this_sz;
  }
  return Rcpp::DataFrame::create(_["estimate_id"]=idx, _["start"]=support, _["beta"]=beta,
                                 _["betahat"]=betahat, _["pseudoweight"]=pseudoweight);
}
       
Rcpp::DataFrame segment_methylation_helper(const Rcpp::DataFrame& df, double tol_val = 0.01);

}

#endif

