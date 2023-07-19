#ifndef BINOMIAL_RATE_DISTRIBUTION_HPP
#define BINOMIAL_RATE_DISTRIBUTION_HPP

/*! \file BinomialRateDistribution.hpp
 * @brief Define BinomialRateDistribution class
 *
 */

#include <Eigen/Core>
#include <Rcpp.h>
#include <cmath> // for std::log

#include "macros.hpp"
#include "Distribution.hpp"
#include "util.hpp" //log(beta) function
#include "Observations.hpp"

namespace MethyLasso {
namespace distribution {

/// @brief Binomial distribution, parametrized in terms of rate
/*! The binomial distribution for a fraction \f$d\f$ of successes,
 * \f$N\f$ trials and success rate \f$\mu\f$ has pmf
 * \f[
 *  p(K|N,\mu) = {N \choose K} \mu^K (1-\mu)^{N-K}
 *  \f]
 * where \f$K=dN\f$. The data to be provided is \f$d\f$ and \f$N\f$ hence the
 * parametrization. The quantity of interest is the rate \f$\mu\f$ such that
 *  \f$E(d)\equiv \mu\f$, or with \f$K\f$ the number of successes, \f$E(K)\equiv N\mu\f$
 *  IRLS does not require \f$N\f$ to be an integer, so we relax that constraint.
 *  
 * - The variance of \f$d\f$ is \f$\mu(1-\mu)/N\f$
 * - In IRLS, the initial guess for \f$\mu\f$ is \f$d\f$
 * - There are no hyperparameters
 */
class BinomialRateDistribution : public Distribution<BinomialRateDistribution> {
public:
  
  BinomialRateDistribution(double unused=0) {}
  
  predictor::MeanVector get_mean_initial_guess_impl(const data::Observations& obs) const {
    return predictor::MeanVector(obs.get_observed());
  }
  
  Eigen::ArrayXd get_variance_impl(const Eigen::ArrayXd& N,
                                   const Eigen::ArrayXd& mu) const {
    return mu*(1-mu)/N;
  }
  
  double get_minus_log_pdf_impl(const Eigen::ArrayXd& d,
                                const Eigen::ArrayXd& N,
                                const Eigen::ArrayXd& m) const {
    Eigen::ArrayXd a = (N+1).log() + util::log_beta(N*(1-d)+1, N*d+1);
    Eigen::ArrayXd b = -(N*d)*(m.log());
    Eigen::ArrayXd c = -(N*(1-d))*((1-m).log());
    /*double ret=0;
    for (unsigned i=0; i<d.size(); ++i) {
       ret += a(i)+b(i)+c(i);
       Rcpp::Rcout << " get_minus_log_pdf d= " << d(i) << " N= " << N(i) << " mu= " << m(i) <<
          " a= " << a(i) << " b= " << b(i) << " c= " << c(i) << " ret= " << ret <<  "\n";
    }*/
    return (a+b+c).sum();
  }
  
  Eigen::ArrayXd get_hyperparams_impl() const { return Eigen::ArrayXd::Constant(1,-1); }
  
  void update_hyperparams_impl(const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&) {}
   void damp_hyperparams_impl(double) {}
  
};


}
}

#endif

