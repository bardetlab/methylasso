#ifndef ESTIMATOR_GROUP_HPP
#define ESTIMATOR_GROUP_HPP

/*! \file EstimatorGroup.hpp
  * @brief Defines the EstimatorGroup class
  */

#include <numeric>

#include "traits.hpp"
#include "DataDesign.hpp"
#include "Estimator.hpp"
#include "macros.hpp"


namespace MethyLasso {
namespace estimator {

/*! @brief The EstimatorGroup estimates a (larger) portion of the Posterior using IRLS
* 
* An EstimatorGroup holds one or several Estimator classes, which have the same
* Leg and Method, but differ in the value of their design
* 
* @tparam Leg the Leg \f$i\f$ of this estimator
* @tparam Method the prior placed on the \f$X_{i\cdot}\f$
* 
*/
template<typename Leg, typename Method>
class EstimatorGroup {
public:
  
  /// @brief Constructs the EstimatorGroup given Data, fixed Config parameters and a DataDesign matrix
  template<typename Data>
  EstimatorGroup(const Data& data, const Config<Leg,Method>& conf, const data::DataDesign& design)
    : estimators_(make_estimators(data, conf, design)), design_(design) {}
  
  /// @brief Update the IRLS mean and weight of each estimator
  void update_binned_residuals(const predictor::IRLSResiduals& z) {
    for (auto&& est : estimators_) est.update_binned_residuals(z);
  }
  
  /// @brief Update the parameters of each estimator
  void update_params() {
    for (auto&& est : estimators_) est.update_params();
  }
  
  /// @brief Damp the parameters of each estimator
  void damp_params(double new_param_weight = 1) {
    for (auto&& est : estimators_) est.damp_params(new_param_weight);
  }
  
  /// @brief Update the hyperparameters of each estimator
  void update_hyperparams() {
    for (auto&& est : estimators_) est.update_hyperparams();
  }
  
  /// @brief Damp the hyperparameters of each estimator
  void damp_hyperparams(double new_param_weight = 1) {
    for (auto&& est : estimators_) est.damp_hyperparams(new_param_weight);
  }
  
  /// @brief Turn on / off the optimization of hyperparameters of each estimator
  void switch_hyperparams_update(bool on = true) {
    for (auto&& est : estimators_) est.switch_hyperparams_update(on);
  }
  
  /// @brief Return the current value of the hyperparameters
  Eigen::ArrayXd get_hyperparams() const {
    unsigned size=0;
    std::vector<Eigen::ArrayXd> tmp;
    for (const auto& est : estimators_) {
      tmp.push_back(est.get_hyperparams());
      size += tmp.back().size();
    }
    Eigen::ArrayXd ret(size);
    unsigned cumsz=0;
    for (const auto& h : tmp) {
      ret.segment(cumsz,h.rows()) = h;
      cumsz += h.rows();
    }
    return ret;
  }
  
  /// @brief Get estimate along support in original data
  predictor::LinearResponse get_linear_response() const { 
    std::vector<Eigen::ArrayXd> tmp;
    for (const auto& est : estimators_) tmp.push_back(est.get_linear_response().get_eta());
    Eigen::ArrayXd ret(Eigen::ArrayXd::Zero(tmp.back().size()));
    for (const auto& h : tmp) ret += h;
    return predictor::LinearResponse(ret);
  }
  
  /// @brief return the sum of minus log prior and hyperprior terms
  double get_minus_log_prior() const {
    double val=0;
    for (const auto& est : estimators_) val += est.get_minus_log_prior();
    return val;
  }
 
  /// @brief return the sum of IRLS target and prior. Lower is better.
  double get_irls_target_with_prior() const {
    double val=0;
    for (const auto& est : estimators_) val += est.get_irls_target_with_prior();
    return val;
  }
  
private:
  template<typename Data>
  std::vector<estimator::Estimator<Leg,Method> > make_estimators(const Data& data,
                                                                 const Config<Leg,Method>& conf,
                                                                 const data::DataDesign& design) const {
    unsigned ncols = design.get_num_estimators();
    std::vector<estimator::Estimator<Leg,Method> > ret;
    ret.reserve(ncols);
    //build Binner for first Estimator
    ret.emplace_back(data, conf, data::DataDesign(design.get_col(0)));
    //remaining Estimators just copy the same binner pointer
    for (unsigned i=1; i<ncols; ++i) {
      ret.emplace_back(data, conf, data::DataDesign(design.get_col(i)), ret[0]);
    }
    return ret;
  }
  
  typedef estimator::Estimator<Leg, Method> estimator_t;
  METHLASSO_GET_SET_OWN(std::vector<estimator_t>, estimators)
  METHLASSO_GET_CONST_OWN(data::DataDesign, design)
    
};

}
}

#endif

