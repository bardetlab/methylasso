#ifndef POSTERIOR_HPP
#define POSTERIOR_HPP

/*! \file Posterior.hpp
 * @brief definition of the Posterior class
 */

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>
#include <utility>

#include "IRLSResiduals.hpp"
#include "macros.hpp"
#include "posterior_policies.hpp"
#include "predictor_helpers.hpp"
#include "Observations.hpp"

namespace MethyLasso {

/*! @brief The Posterior is the top-level class to perform IRLS estimation
 *
 * For a bit of context, see \ref index "the main page"
 * 
 * @tparam Distribution the distribution of the likelihood terms (must be all the same)
 * @tparam Link the link function used
 * @tparam IRLSPolicy how the main parameters are updated
 * @tparam HyperPolicy how the distribution and prior hyperparameters are updated
 * 
 */
template<typename Distribution, typename Link, typename IRLSPolicy = BasicIRLSPolicy,
         typename PriorPolicy = BasicPriorPolicy, typename HyperPolicy = BasicHyperPolicy>
class Posterior : private IRLSPolicy, private PriorPolicy, private HyperPolicy {
public:
  
  /// @brief Constructs the Posterior for a given distribution and link function
  /// @param tol_val is the tolerance on \f$\mu\f$ at which convergence is declared 
  /// @param maxiter is the maximum number of iterations performed
  /// @param t remaining parameters forwarded to the Distribution constructor
  template<class... T>
  Posterior(double tol_val, unsigned maxiter, T&&... t)
    : distribution_(std::forward<T>(t)...), tol_val_(tol_val), maxiter_(maxiter),
      converged_(false) {
    //Rcpp::Rcout << "Init Posterior with tol_val=" << tol_val << " and maxiter=" << maxiter << "\n";
  }

  /// @brief Initialize IRLS summaries wihin each estimator to a sensible initial value
  /// Sets zhat and weight in the Summarizer of each Estimator
  /// Sets the converged flag to false
  template<class... Est>
  void initial_guess(const data::Observations& data, Est&&... est);
  
  /// @brief perform IRLS iterations on the estimators provided
  template<class... Est>
  predictor::MeanVector iterate(const data::Observations& data, Est&&... est);
  
  /// @brief compute and return the minus log of posterior, likelihoods and priors (in order)
  template<class... Est>
  Eigen::Array3d get_minus_log_posterior(const data::Observations& data, Est&&... est) const;
  
  /// @brief return the current estimate of the posterior mean \f$\mu = g^{-1}(\phi)\f$
  template<class... T>
  predictor::MeanVector get_mu(T&&... t) const;
  
  /// @brief return distribution hyperparameters
  Eigen::ArrayXd get_distribution_hyperparams() const { return distribution_.get_hyperparams(); }
  
  /// @brief return prior hyperparameters
  template<class Est, class... T> Eigen::ArrayXd get_prior_hyperparams(const Est& est, T&&...) const;
  template<class Est> Eigen::ArrayXd get_prior_hyperparams(const Est& est) const;
  
  /// @brief update all parameters
  /// If successful, will set the converged flag to true
  template<class... Est>
  predictor::MeanVector maximize_posterior(unsigned nouter, bool optimize_prior,
                                    bool optimize_hyper, const data::Observations& data, Est&&... est);
  
  /// @brief whether the maximum posterior has been reached or not
  bool has_converged() const { return get_converged(); }
  
private:

  //parameter update: call update_params() for each estimator
  template<class Est, class... T> void update_all_params(Est&, T&&...);
  template<class Est> void update_all_params(Est&);
  
  //parameter update: apply damping policy
  template<class Data, class... Est>
  void adapt_all_params(unsigned step, double old_mlogp, const Data& data, Est&&... est);
  
  //parameter update: call damp_params() for each estimator
  template<class Est, class... T> void damp_all_params(Est&, T&&...);
  template<class Est> void damp_all_params(Est&);
  
  /// @brief update prior hyperparameters
  template<class Est, class... T> void update_prior_hyperparams(Est&, T&&...);
  template<class Est> void update_prior_hyperparams(Est&);
  
  //parameter update: apply damping policy
  template<class Data, class... Est>
  void adapt_prior_hyperparams(unsigned step, double old_mlogp, const Data& data, Est&&... est);
  
  /// @brief damp prior hyperparameters
  template<class Est, class... T> void damp_prior_hyperparams(Est&, T&&...);
  template<class Est> void damp_prior_hyperparams(Est&);
  
  //hyperparameter update: apply damping policy
  template<class Data, class... Est>
  void adapt_distribution_hyperparams(unsigned step, double old_mlogp, const Data& data, Est&&... est);
  
  /// @brief damp distribution hyperparameters
  void damp_distribution_hyperparams();
  
  //binned_residuals update: call update_summaries(z) for each estimator
  template<class Est, class... T> void update_all_binned_residuals(const predictor::IRLSResiduals&, Est&, T&&...);
  template<class Est> void update_all_binned_residuals(const predictor::IRLSResiduals&, Est&);
  
  //get_eta: return sum of all get_linear_response() calls
  template<class Est, class... T> predictor::LinearResponse get_eta(const Est&, T&&...) const;
  template<class Est> predictor::LinearResponse get_eta(const Est&) const;
  
  //get_minus_log_prior: call get_minus_log_prior() for each estimator
  template<class Est, class... T> double get_minus_log_prior(const Est&, T&&...) const;
  template<class Est> double get_minus_log_prior(const Est&) const;
  
  //get_irls_target_with_prior: call get_irls_target_with_prior() for each estimator
  template<class Est, class... T> double get_irls_target_with_prior(const Est&, T&&...) const;
  template<class Est> double get_irls_target_with_prior(const Est&) const;
  
private:
  Distribution distribution_;
  
  METHLASSO_GET_SET_SIMPLE(double, tol_val)
  METHLASSO_GET_SET_SIMPLE(unsigned, maxiter)
  METHLASSO_GET_SET_SIMPLE(bool, converged)
  
};

#include "Posterior.ipp"

}

#endif

