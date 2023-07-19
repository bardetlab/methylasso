#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP

/*! \file Estimator.hpp
  * @brief Defines the Estimator class
  */

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>
#include <algorithm> //std::min
#include <memory> //shared_ptr

#include "IRLSResiduals.hpp"
#include "macros.hpp"
#include "predictor_helpers.hpp"
#include "binned_residuals_helpers.hpp"
#include "traits.hpp"
#include "Design.hpp"
#include "DataDesign.hpp"

#include "Estimate.hpp"
#include "estimate_helpers.hpp"

#include "Pseudodata.hpp"
#include "pseudodata_helpers.hpp"
#include "params_helpers.hpp"

namespace MethyLasso {
namespace estimator {

/*! @brief The Estimator estimates a portion of the Posterior using IRLS
* 
* @tparam Leg the Leg \f$i\f$ of this estimator
* @tparam Method the prior placed on the \f$X_{i\cdot}\f$
* 
*/
template<typename Leg, typename Method>
class Estimator {
public:
  
  /// @brief Constructs the Estimator given Data and fixed Config parameters
  template<typename Data>
  Estimator(const Data& data, const Config<Leg,Method>& conf, const data::DataDesign& data_design)
    : binner_ptr_(make_binner_ptr<Leg,Method>(data,conf)), binned_residuals_(),
      design_(make_design(get_binner(), conf, data_design)),
      pseudodata_(), params_(design_,conf), old_params_(params_), hyper_(conf), old_hyper_(hyper_),
      fitter_(design_,conf), update_hyper_(true), data_design_(data_design) {
    //check for user interrupt
    Rcpp::checkUserInterrupt();
  }
  
  /// @brief Constructs the Estimator given Data and fixed Config parameters, reusing another Estimator's Binner
  /*! @note you should not use this constructor unless you know what you're doing. Use EstimatorGroup instead. */
  template<typename Data>
  Estimator(const Data& data, const Config<Leg,Method>& conf, const data::DataDesign& data_design, const Estimator<Leg,Method>& other)
    : binner_ptr_(other.get_binner_ptr()), binned_residuals_(),
      design_(make_design(get_binner(), conf, data_design)),
      pseudodata_(), params_(design_,conf), old_params_(params_), hyper_(conf), old_hyper_(hyper_),
      fitter_(design_,conf), update_hyper_(true), data_design_(data_design) {
          //check for user interrupt
          Rcpp::checkUserInterrupt();
      }
  
  /// @brief Update the IRLS mean and weight
  void update_binned_residuals(const predictor::IRLSResiduals& z);
  
  /// @brief Update the parameters using the IRLS summaries
  void update_params();
  
  /// @brief Perform damping of parameters, mixing old and new parameter estimates
  void damp_params(double new_param_weight = 1);
  
  /// @brief Update the hyperparameters using the current estimates of the parameters
  void update_hyperparams();
  
  /// @brief Perform damping of hyperparameters
  void damp_hyperparams(double new_param_weight = 1);
  
  /// @brief Turn on / off the optimization of hyperparameters. By default, is on.
  void switch_hyperparams_update(bool on = true) { set_update_hyper(on); }
  
  /// @brief Return the current value of the hyperparameters
  Eigen::ArrayXd get_hyperparams() const { return get_hyper().get_all(); }
  
  /// @brief Get estimate along support in binned data
  binned::Estimate get_estimate() const { 
    binned::Estimate est = binned::get_estimate(get_params(), get_design());
    //print debug info
    if (estimator::EstimatorTraits<Leg,Method>::debug_response) {
      Rcpp::Rcout << "   estimator " << (this) << " get_estimate(): "
                  << est.get_phi().tail(std::min<int>(10,est.get_phi().rows())).transpose() << "\n";
    }
    return est;
  }
  
  /// @brief Get estimate along support in original data
  predictor::LinearResponse get_linear_response() const {
    auto lr = predictor::get_linear_response(get_estimate(), get_binner());
    //print debug info
    if (estimator::EstimatorTraits<Leg,Method>::debug_response) {
      Rcpp::Rcout << "  estimator " << (this) << " get_linear_response(): "
                  << lr.get_eta().tail(std::min<int>(10,lr.get_eta().rows())).transpose() << "\n";
    }
    return lr;
  }
  
  /// @brief return the sum of minus log prior and hyperprior terms
  double get_minus_log_prior() const {
    return get_fitter().get_minus_log_prior(get_params(), get_hyper());
  }
  
  /// @brief return the sum of IRLS target and prior. Lower is better.
  double get_irls_target_with_prior() const {
    double prior = get_minus_log_prior();
    Eigen::ArrayXd vec = get_pseudodata().get_betahat() - get_params().get_beta();
    double target = 0.5 * (vec * get_pseudodata().get_pseudoweight() * vec).sum();
    return target + prior;
  }
  
public:
  typedef Fitter<Method, typename EstimatorTraits<Leg,Method>::library > fitter_t;
  
  const binned::Binner& get_binner() const { return *get_binner_ptr(); }
  
  METHLASSO_GET_CONST_OWN(std::shared_ptr<binned::Binner>, binner_ptr) //const and agnostic to Leg and Method
  METHLASSO_GET_SET_OWN(binned::BinnedResiduals, binned_residuals) //agnostic to Leg and Method
    
  METHLASSO_GET_CONST_OWN(params::Design, design)
  METHLASSO_GET_SET_OWN(params::Pseudodata, pseudodata)
  METHLASSO_GET_SET_OWN(params::Params<Method>, params)
  METHLASSO_GET_SET_OWN(params::Params<Method>, old_params)
  METHLASSO_GET_SET_OWN(params::Hyperparams<Method>, hyper)
  METHLASSO_GET_SET_OWN(params::Hyperparams<Method>, old_hyper)
    
  METHLASSO_GET_CONST_OWN(fitter_t, fitter)
  METHLASSO_GET_SET_SIMPLE(bool, update_hyper)
  
  METHLASSO_GET_CONST_OWN(data::DataDesign, data_design)
    
};

#include "Estimator.ipp"

}
}

#endif

