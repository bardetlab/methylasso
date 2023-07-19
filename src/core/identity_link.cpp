#include "link_functions.hpp"

namespace MethyLasso {
namespace link {

template<>
predictor::LinearResponse get_eta<Identity>(const predictor::MeanVector& mu) {
  return predictor::LinearResponse(mu.get_mu());
}

template<>
predictor::MeanVector get_mu<Identity>(const predictor::LinearResponse& eta) {
  return predictor::MeanVector(eta.get_eta());
}

template<>
predictor::MeanVector restrict_to_support<Identity>(const predictor::MeanVector& mu0) { return mu0; }

template<>
predictor::IRLSResiduals transform_initial_guess<Identity>(const predictor::Moments& mo) {
  Eigen::ArrayXd residuals(mo.get_mean()); // assume initial guess of parameters is zero
  Eigen::ArrayXd weights(1/mo.get_variance());
  /*Rcpp::Rcout << "initial guess:\n";
   Rcpp::Rcout << "   mu0 " << mo.get_mean().head(3).transpose() << "\n";
   Rcpp::Rcout << "   V0 " << mo.get_variance().head(3).transpose() << "\n";
   Rcpp::Rcout << "   residuals " << residuals.head(3).transpose() << "\n";
   Rcpp::Rcout << "   weights " << weights.head(3).transpose() << "\n";*/
  return predictor::IRLSResiduals{residuals,weights};
}

template<>
predictor::IRLSResiduals get_irls_residuals<Identity>(const data::Observations& obs, const predictor::Moments& mo) {
  Eigen::ArrayXd residuals(obs.get_observed()-mo.get_mean()); // return z
  Eigen::ArrayXd weights(1/mo.get_variance());
  return predictor::IRLSResiduals{residuals,weights};
}

}
}


