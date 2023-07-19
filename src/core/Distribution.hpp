#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

/*! \file Distribution.hpp
 * @brief Define base Distribution class
 *
 */

#include <Eigen/Core>
#include <Rcpp.h>

#include "MeanVector.hpp"

namespace MethyLasso {
namespace distribution {

/// @brief Class which describes the interface to a disrtribution
template<class Derived>
class Distribution {
public:
  
  /// @brief Compute a pointwise estimate of the mean, used as initial guess for IRLS
  /* Ensures that the value returned is within the support of the distribution
   */
  template<class... T>
  predictor::MeanVector get_mean_initial_guess(T&&... t) const {
    return static_cast<const Derived*>(this)->get_mean_initial_guess_impl(std::forward<T>(t)...);
  }
  
  /// @brief Compute the value of the variance of the distribution, given the mean
  /* Make sure the variance is defined and >0 at that value of the mean, otherwise behaviour is undefined
   */
  template<class... T>
  Eigen::ArrayXd get_variance(T&&... t) const {
    return static_cast<const Derived*>(this)->get_variance_impl(std::forward<T>(t)...);
  }
  
  /// @brief Compute the sum of minus log pdf on each data point
  template<class... T>
  double get_minus_log_pdf(T&&... t) const {
    return static_cast<const Derived*>(this)->get_minus_log_pdf_impl(std::forward<T>(t)...);
  }
  
  /// @brief Update any parameters of that distribution which are held fixed for IRLS
  /*! For this step, we assume a flat prior for the hyperparameter */
  template<class... T>
  void update_hyperparams(T&&... t) {
    return static_cast<Derived*>(this)->update_hyperparams_impl(std::forward<T>(t)...);
  }
  
  /// @brief Damp any parameters of that distribution which are held fixed for IRLS
  /*! For this step, we assume a flat prior for the hyperparameter */
  template<class... T>
  void damp_hyperparams(double new_param_weight, T&&... t) {
    return static_cast<Derived*>(this)->damp_hyperparams_impl(new_param_weight,std::forward<T>(t)...);
  }
  
  Eigen::ArrayXd get_hyperparams() const {
    return static_cast<const Derived*>(this)->get_hyperparams_impl();
  }
  
  
};


}
}

#endif

