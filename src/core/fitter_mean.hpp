#ifndef FITTER_MEAN_HPP
#define FITTER_MEAN_HPP

/*! \file fitter_mean.hpp
 * @brief Defines Params, Hyperparams and Fitter for the Mean method
 */

#include <RcppEigen.h>
using namespace Rcpp;
#include <Eigen/Core>

#include "macros.hpp"
#include "traits.hpp"

#include "Design.hpp"

namespace MethyLasso {
namespace binned {
class Binner;
}
namespace params {

// class that holds parameters that were output by a Mean fit
template<>
class Params<Mean> {
public:
  template<typename Config>
  Params(const params::Design& design, const Config&) : beta_(Eigen::ArrayXd::Zero(design.get_nparams())) {}
  Params(const Eigen::ArrayXd& beta) : beta_(beta) {}
  
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, beta)
    
};

// class that holds hyperparameters that were output by a Mean fit
template<>
class Hyperparams<Mean> {
public:
  Hyperparams() {}
  template<typename Config>
  Hyperparams(const Config&) : Hyperparams() {}
  Eigen::ArrayXd get_all() const { return Eigen::ArrayXd::Constant(0,0); }
  
};

}
namespace estimator {

template<typename Library>
class Fitter<Mean,Library> {
public:
  template<typename Config>
  Fitter(const params::Design&, const Config&) {}
   
  //update beta given betahat and pseudoweight
  params::Params<Mean> fit_params(const params::Pseudodata& pseudo, const params::Hyperparams<Mean>& hyper) const {
      //simply set beta = betahat
      Eigen::ArrayXd betahat = pseudo.get_betahat(); //essential, because the std::vector uses that memory location
      Eigen::ArrayXd pseudoweight = pseudo.get_pseudoweight();
      return params::Params<Mean>{betahat};
  }

  params::Params<Mean> damp_params(const params::Params<Mean>& new_p,
                                    const params::Params<Mean>& old_p,
                                    double new_param_weight) const {
    auto damped_p = (old_p.get_beta().size()>0)
                    ? params::Params<Mean>(new_param_weight*new_p.get_beta() + (1-new_param_weight)*old_p.get_beta())
                    : new_p;
    return damped_p;
  }
  
  params::Hyperparams<Mean> fit_hyperparams(const params::Params<Mean>& par) const {
      return params::Hyperparams<Mean>{};
  }
  
  params::Hyperparams<Mean> damp_hyperparams(const params::Hyperparams<Mean>&,
                                              const params::Hyperparams<Mean>&, double) const {
      return params::Hyperparams<Mean>{};
  }
  
  /// @brief return the sum of minus log prior and hyperprior terms
  double get_minus_log_prior(const params::Params<Mean>&,
                             const params::Hyperparams<Mean>&) const { return 0; }
  
};

}
}

#endif

