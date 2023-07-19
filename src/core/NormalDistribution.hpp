#ifndef NORMAL_DISTRIBUTION_HPP
#define NORMAL_DISTRIBUTION_HPP

/*! \file NormalDistribution.hpp
 * @brief Define NormalDistribution class
 */

#include <Eigen/Core>
#include <Rcpp.h>
#include <cmath> // for M_PI and std::sqrt

#include "macros.hpp"
#include "Distribution.hpp"
#include "Observations.hpp"

namespace MethyLasso {
namespace distribution {

/// @brief Normal distribution
/*! The normal distribution for observation \f$y\f$, weight \f$w\f$,
 *  mean \f$\mu\f$ and standard deviation parameter \f$\sigma\f$ has pdf
 * \f[
 * p(y|\mu,w,\sigma) = \frac{1}{(\sqrt{2\pi}\sigma)^w}\exp\left(-w\frac{(y-\mu)^2}{2\sigma^2}\right)
 * \f]
 * 
 * - The variance is \f$\sigma^2/w\f$
 * - In IRLS, the initial guess for \f$\mu\f$ is \f$y\f$
 * - \f$\sigma\f$ is considered a hyperparameter
 * 
 */
class NormalDistribution : public Distribution<NormalDistribution> {
public:
  NormalDistribution(double sigma = 1) : sigma_(sigma), old_sigma_(sigma_) {
    //Rcpp::Rcout << "Init normal distribution with sigma=" << sigma_ << '\n';
  }
  
  predictor::MeanVector get_mean_initial_guess_impl(const data::Observations& obs) const {
    return predictor::MeanVector(obs.get_observed());
  }
  
  Eigen::ArrayXd get_variance_impl(const Eigen::ArrayXd& wt,
                                   const Eigen::ArrayXd&) const {
    return SQUARE(get_sigma())/wt;
  }
  
  double get_minus_log_pdf_impl(const Eigen::ArrayXd& y,
                                const Eigen::ArrayXd& wt,
                                const Eigen::ArrayXd& mu) const {
    double sigma = get_sigma();
    return ( wt/2.*( ((y-mu)/sigma).square() + std::log(2*M_PI*sigma*sigma) ) ).sum();
  }
  
  void update_hyperparams_impl(const Eigen::ArrayXd& y,
                               const Eigen::ArrayXd& wt,
                               const Eigen::ArrayXd& mu) {
    set_old_sigma(get_sigma());
    double sigma = std::sqrt( ((y-mu).square()*wt).sum()/wt.sum() );
    set_sigma(sigma);
    //Rcpp::Rcout << "   update_hyperparams(): new sigma = " << get_sigma() << "\n";
  }

  void damp_hyperparams_impl(double new_param_weight) {
    set_sigma(new_param_weight*get_sigma()+(1-new_param_weight)*get_old_sigma());
    //Rcpp::Rcout << "   update_hyperparams(): new sigma = " << get_sigma() << "\n";
  }
  
  Eigen::ArrayXd get_hyperparams_impl() const { return Eigen::ArrayXd::Constant(1,get_sigma()); }
  
  METHLASSO_GET_SET_SIMPLE(double, sigma)
  METHLASSO_GET_SET_SIMPLE(double, old_sigma)
};


}
}

#endif

