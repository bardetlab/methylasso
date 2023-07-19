#include "link_functions.hpp"

namespace MethyLasso {
namespace link {

template<>
predictor::LinearResponse get_eta<Logit>(const predictor::MeanVector& mv) {
  auto mu = mv.get_mu();
  return predictor::LinearResponse((mu/(1-mu)).log()); 
}

template<>
predictor::MeanVector get_mu<Logit>(const predictor::LinearResponse& lr) {
  auto eta = lr.get_eta();
  return predictor::MeanVector(1/( (-eta).exp() + 1));
}

template<>
predictor::MeanVector restrict_to_support<Logit>(const predictor::MeanVector& mu0) {
  return predictor::MeanVector(mu0.get_mu().cwiseMin(1-1e-10).cwiseMax(1e-10));
}

template<>
predictor::IRLSResiduals transform_initial_guess<Logit>(const predictor::Moments& mo) {
  auto mu0 = mo.get_mean();
  auto V0 = mo.get_variance();
  Eigen::ArrayXd residuals(get_eta<Logit>(mu0).get_eta()); // assume initial guess of parameters is zero
  Eigen::ArrayXd weights( ( mu0*(1-mu0) ).square()/V0 ) ;
  return predictor::IRLSResiduals{residuals,weights};
}

template<>
predictor::IRLSResiduals get_irls_residuals<Logit>(const data::Observations& obs, const predictor::Moments& mo) {
  auto mu = mo.get_mean();
  auto Vk = mo.get_variance();
  Eigen::ArrayXd residuals( (obs.get_observed()-mu)/(mu*(1-mu)) ); // return z
  Eigen::ArrayXd weights( ( mu*(1-mu) ).square()/Vk );
  return predictor::IRLSResiduals{residuals,weights};
}


}
}


