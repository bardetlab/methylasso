#ifndef BETA_BINOMIAL_RATE_DISTRIBUTION_HPP
#define BETA_BINOMIAL_RATE_DISTRIBUTION_HPP

/*! \file BetaBinomialRate.hpp
  * @brief Define BetaBinomialRate class
  *
  */

#include <Eigen/Core>
#include <Rcpp.h>
#include <cmath> // for std::log
#include <boost/math/tools/minima.hpp> // to find hyperparameter
#include <boost/cstdint.hpp>

#include "macros.hpp"
#include "Distribution.hpp"
#include "util.hpp" //log(beta) function
#include "Observations.hpp"

namespace MethyLasso {
namespace distribution {

/// @brief Beta-binomial distribution, parametrized in terms of rate
/*! The binomial distribution for \f$N\f$ trials, \f$dN\f$ successes
* and success rate \f$p\f$ has pmf
* \f[
*  p(dN|N,p) = {N \choose dN} p^{dN} (1-p)^{N(1-d)}
*  \f]
*  Under this distribution, \f$E(d)=p\f$ and \f$V(d)=p(1-p)/N\f$.
*  The beta distribution of \f$p\f$ with mean \f$\mu\f$ and sample size \f$\nu\f$ has pdf
*  \f[
*  p(p|\mu,\nu) = \frac{p^{\nu\mu-1}(1-p)^{\nu(1-\mu)-1}}{B(\nu\mu,\nu(1-\mu))}
*  \f]
*  where \f$B\f$ is the beta function.
*  The beta-binomial distribution is obtained by placing a beta prior on \f$p\f$ when \f$K\f$
*  follows a binomial distribution. Marginalizing \f$p\f$ yields the pdf
*  \f[
*  p(dN|N,\mu,\nu) = {N \choose dN} \frac{B(dN+\nu\mu,N(1-d)+\nu(1-\mu))}{B(\nu\mu,\nu(1-\mu))}
*  \f]
*  
* The data to be provided is \f$d\f$ and \f$N\f$. The quantity of interest is the rate \f$\mu\f$.
*  
*  - \f$E(d)\equiv \mu\f$
*  - \f$V(d)=\frac{\mu(1-\mu)}{N}\frac{N+\nu}{1+\nu}\f$. The distribution has overdispersion, and
*    reverts to the binomial distribution when \f$\nu\to+\infty\f$
*  - In IRLS, the initial guess for \f$\mu\f$ is \f$d\f$
*  - The parameter \f$\nu\f$ is treated as a hyperparameter. It is found by finding the root of the derivative of
*    \f$-\log p(dN|N,\mu,\nu)\f$ :
*    \f[ \frac{D-\log p}{D\nu} = \mu\left(\psi(\mu\nu)-\psi(dN+\mu\nu)\right)
*                              + (1-\mu)\left(\psi(\nu(1-\mu))-\psi(N(1-d)+\nu(1-\mu))\right)
*                              - \psi(\nu) + \psi(N+\nu)
*    \f]
*    where \f$\psi\f$ is the digamma function.
*/
class BetaBinomialRateDistribution : public Distribution<BetaBinomialRateDistribution> {
public:
  
  BetaBinomialRateDistribution(double nu = 1) : nu_(nu), old_nu_(nu_) {}
  
  predictor::MeanVector get_mean_initial_guess_impl(const data::Observations& obs) const {
    return predictor::MeanVector(obs.get_observed());
  }
  
  Eigen::ArrayXd get_variance_impl(const Eigen::ArrayXd& N,
                                   const Eigen::ArrayXd& mu) const {
    return mu*(1-mu)/N * (N+get_nu())/(1+get_nu());
  }
  
  double get_minus_log_pdf_impl(const Eigen::ArrayXd& d,
                                const Eigen::ArrayXd& N,
                                const Eigen::ArrayXd& m) const {
    double nu=get_nu();
    Eigen::ArrayXd a = (N+1).log() + util::log_beta(N*(1-d)+1, N*d+1);
    Eigen::ArrayXd b = -util::log_beta(N*d+m*nu, N*(1-d)+nu*(1-m));
    Eigen::ArrayXd c = util::log_beta(m*nu, nu*(1-m)); 
    /*double ret=0;
    for (unsigned i=0; i<d.size(); ++i) {
    ret += a(i)+b(i)+c(i);
    Rcpp::Rcout << " get_minus_log_pdf d= " << d(i) << " N= " << N(i) << " mu= " << m(i) <<
      " a= " << a(i) << " b= " << b(i) << " c= " << c(i) << " ret= " << ret <<  "\n";
    }*/
    for (unsigned i=0; i<a.size(); ++i)
      if (std::isnan(a(i)+b(i)+c(i))) {
      Rcpp::Rcout << "NaN found at i= " << i << " d= " << d(i) << " N= " << N(i) << " m= " << m(i)
                  << " a= " << a(i) << " b= " << b(i) << " c= " << c(i) << " nu= " << nu << " \n";
      Rcpp::Rcout << "   m*nu= " << m(i)*nu << " (1-m)*nu= " << (1-m(i))*nu << "\n";
      break;
      }
    return (a+b+c).sum();
  }
   
  //adapted from boost doc for brent_minima
  void update_hyperparams_impl(const Eigen::ArrayXd& d,
                               const Eigen::ArrayXd& N,
                               const Eigen::ArrayXd& m) {
    auto target = [&] (double nu) {
      Eigen::ArrayXd b = -util::log_beta(N*d+m*nu, N*(1-d)+nu*(1-m));
      Eigen::ArrayXd c = util::log_beta(m*nu, nu*(1-m)); 
      return (b+c).sum();
    };
    int bits = std::numeric_limits<double>::digits;
    boost::uintmax_t maxiter = 100;
    double bracket = 1e6; //look around current nu in this range in log scale
    double nu = get_nu();
    set_old_nu(nu);
    std::pair<double, double> r = boost::math::tools::brent_find_minima(target, nu/bracket, nu*bracket, bits, maxiter);
    //Rcpp::Rcout << "   update_hyperparams: nu at minimum = " << r.first << ", -log p(" << r.first << ") = " << r.second << std::endl;
    //Rcpp::Rcout << "   update_hyperparams(): use damp= " << damp << " and factor= " << new_param_weight << "\n";
    set_nu(r.first);
    if (get_nu() <= nu/bracket || get_nu() >= nu*bracket) Rcpp::Rcout << "Warning: nu optimization hit boundary\n";
    //Rcpp::Rcout << "   update_hyperparams(): new nu = " << get_nu() << "\n";
  }
   
  void damp_hyperparams_impl(double new_param_weight) {
      set_nu(new_param_weight*get_nu()+(1-new_param_weight)*get_old_nu());
    //Rcpp::Rcout << "   update_hyperparams(): new nu = " << get_nu() << "\n";
  }
  
  Eigen::ArrayXd get_hyperparams_impl() const { return Eigen::ArrayXd::Constant(1,get_nu()); }
  
  METHLASSO_GET_SET_SIMPLE(double, nu)
  METHLASSO_GET_SET_SIMPLE(double, old_nu)
  
};

}
}

#endif

