#include "binned_residuals_helpers.hpp"
#include <Rcpp.h>

#include <Eigen/Core>

namespace MethyLasso {
namespace binned {

BinnedResiduals bin_irls_residuals(const predictor::IRLSResiduals& z, const binned::Binner& binner) {
  //compute weight
  Eigen::ArrayXd weight_sum = binner.bin(z.get_weights()); //equivalent of diagonal of B_l W^k B_l^\top
  //compute zhat
  Eigen::ArrayXd zhat = ( binner.bin(z.get_residuals() * z.get_weights()) / weight_sum );
  //fix bins with weight_sum == 0 (i.e. nobs == 0) by setting aggregate z to 0
  zhat = (weight_sum==0).select(weight_sum,zhat);
  return BinnedResiduals{zhat,weight_sum};
}

}
}

