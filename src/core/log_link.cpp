#include "link_functions.hpp"

namespace MethyLasso {
namespace link {

template<>
predictor::LinearResponse get_eta<Logarithmic>(const predictor::MeanVector& mu) {
  return predictor::LinearResponse(mu.get_mu().log());
}

template<>
predictor::MeanVector get_mu<Logarithmic>(const predictor::LinearResponse& eta) {
  return predictor::MeanVector(eta.get_eta().exp());
}

template<>
predictor::MeanVector restrict_to_support<Logarithmic>(const predictor::MeanVector& mu0) {
  return predictor::MeanVector(mu0.get_mu().cwiseMax(1e-3));
}

template<>
predictor::IRLSResiduals transform_initial_guess<Logarithmic>(const predictor::Moments& mo) {
  Eigen::ArrayXd residuals(get_eta<Logarithmic>(mo.get_mean()).get_eta()); // assume initial guess of parameters is zero
  Eigen::ArrayXd weights( mo.get_mean().square()/mo.get_variance() );
  //Rcpp::Rcout << " link mu0= " << mo.get_mean().transpose().head(3) << "\n";
  //Rcpp::Rcout << " link V0= " << mo.get_variance().transpose().head(3) << "\n";
  return predictor::IRLSResiduals{residuals,weights};
}

template<>
predictor::IRLSResiduals get_irls_residuals<Logarithmic>(const data::Observations& obs, const predictor::Moments& mo) {
  Eigen::ArrayXd residuals( obs.get_observed()/mo.get_mean()-1 ); // return z
  Eigen::ArrayXd weights( mo.get_mean().square()/mo.get_variance() );
  return predictor::IRLSResiduals{residuals,weights};
}

}
}


