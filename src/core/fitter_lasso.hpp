#ifndef FITTER_LASSO_HPP
#define FITTER_LASSO_HPP

/*! \file fitter_lasso.hpp
 * @brief Defines Params, Hyperparams and Fitter for the Lasso method
 */

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "macros.hpp"
#include "traits.hpp"

#include "FusedLassoGaussianEstimator.hpp"
#include "Design.hpp"

namespace MethyLasso {
namespace binned {
class Binner;
}
namespace params {

// class that holds parameters that were output by a Lasso fit
template<>
class Params<Lasso> {
public:
  template<typename Config>
  Params(const params::Design& design, const Config&) : beta_(Eigen::ArrayXd::Zero(design.get_nparams())) {}
  Params(const Eigen::ArrayXd& beta) : beta_(beta) {}
  
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, beta)
    
};

// class that holds hyperparameters that were output by a Lasso fit
template<>
class Hyperparams<Lasso> {
public:
  Hyperparams(double lambda2) : lambda2_(lambda2) {}
  template<typename Config>
  Hyperparams(const Config& conf) : Hyperparams(conf.get_lambda2()) {}
  Eigen::ArrayXd get_all() const { return Eigen::ArrayXd::Constant(1,get_lambda2()); }
  METHLASSO_GET_SET_SIMPLE(double, lambda2)
  
};

}
namespace estimator {

template<typename Library>
class Fitter<Lasso,Library> {
public:
  template<typename Config>
  Fitter(const params::Design& design, const Config& conf) : 
    flo_(design.get_nparams(), conf.get_tol_val(), conf.get_max_lambda2()), lambda1_(conf.get_lambda1()) {}
   
  //update beta given betahat and pseudoweight
  params::Params<Lasso> fit_params(const params::Pseudodata& pseudo, const params::Hyperparams<Lasso>& hyper) const;

  params::Params<Lasso> damp_params(const params::Params<Lasso>& new_p,
                                    const params::Params<Lasso>& old_p,
                                    double new_param_weight) const {
    auto damped_p = (old_p.get_beta().size()>0)
                    ? params::Params<Lasso>(new_param_weight*new_p.get_beta() + (1-new_param_weight)*old_p.get_beta())
                    : new_p;
    return damped_p;
  }
  
  //update lambda2 given beta
  params::Hyperparams<Lasso> fit_hyperparams(const params::Params<Lasso>& par) const;
  
  params::Hyperparams<Lasso> damp_hyperparams(const params::Hyperparams<Lasso>& new_p,
                                              const params::Hyperparams<Lasso>& old_p,
                                              double new_param_weight) const {
    auto damped_p = params::Hyperparams<Lasso>(new_param_weight*new_p.get_lambda2() + (1-new_param_weight)*old_p.get_lambda2());
    return damped_p;
  }
  
  
  /// @brief return the sum of minus log prior and hyperprior terms
  double get_minus_log_prior(const params::Params<Lasso>& par,
                             const params::Hyperparams<Lasso>& hyper) const;
  
private:
  //used to fit parameters. Mutable because internal caching can occur.
  mutable base::FusedLassoGaussianEstimator<Library> flo_;
  METHLASSO_GET_CONST_SIMPLE(double, lambda1)
};

#include "fitter_lasso.ipp"

}
}

#endif

