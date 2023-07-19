#ifndef METHYLATION_DETECTION_HPP
#define METHYLATION_DETECTION_HPP

#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "MethData.hpp"
#include "Posterior.hpp"
#include "offset.hpp"
#include "posterior_policies.hpp"
#include "MeanVector.hpp"
#include "DataDesign.hpp"
#include "EstimatorGroup.hpp"
#include "methylation_helpers.hpp"
#include "data_design_helpers.hpp"
#include "binned_genome.hpp"

namespace MethyLasso {

template<class Distribution, class Link, class Leg = BinnedGenome,
         class IRLSPolicy = BasicIRLSPolicy, class PriorPolicy = BasicPriorPolicy, class HyperPolicy = BasicHyperPolicy>
Rcpp::List methylation_detection(const MethData& obs, const data::DataDesign& design, const estimator::Config<Leg, Lasso>& mconf,
                                 bool optimize_lambda2=true, double hyper=1, bool optimize_hyper=true) {
  //build Posterior
  //Rcpp::Rcout << "build Posterior\n";
  Posterior<Distribution, Link, IRLSPolicy, PriorPolicy, HyperPolicy> post(mconf.get_tol_val(), mconf.get_max_iter(), hyper);
  //
  //build estimators
  //Rcpp::Rcout << "build offset estimator\n";
  OffsetEstimator offset(obs, OffsetConfig(), data::DataDesign());
  //Rcpp::Rcout << "build genomic estimators\n";
  estimator::EstimatorGroup<Leg, Lasso> meth(obs, mconf, design);
  meth.switch_hyperparams_update(optimize_lambda2);
  //
  //Initialize IRLS
  //Rcpp::Rcout << "Initial guess\n";
  post.initial_guess(obs, meth, offset);
  //
  //run IRLS iterations
  //Rcpp::Rcout << "Begin IRLS iterations\n";
  predictor::MeanVector mu = post.maximize_posterior(mconf.get_nouter(), optimize_lambda2,
                                                     optimize_hyper, obs, meth, offset);
  //
  //return transformed data
  auto offset_final = offset.get_params().get_beta();
  auto prior_final = meth.get_hyperparams();
  auto hyper_final = post.get_distribution_hyperparams();
  bool converged = post.has_converged();
  auto est_df = get_estimate_dataframe(meth);
  return Rcpp::List::create(_["mu"]=mu.get_mu(), _["lambda2"]=prior_final, _["hyper"]=hyper_final,
                            _["converged"]=converged, _["offset"]=offset_final, _["estimates"]=est_df);
}

}

#endif

