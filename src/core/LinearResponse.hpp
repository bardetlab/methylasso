#ifndef LINEAR_RESPONSE_HPP
#define LINEAR_RESPONSE_HPP

/*! \file LinearResponse.hpp
 * @brief Defines LinearResponse class
 */

#include <Eigen/Core>

#include "macros.hpp"

namespace MethyLasso {
namespace predictor {

//! @brief LinearResponse holds \f$\eta = B^\top X \beta \f$
struct LinearResponse {
  LinearResponse(const Eigen::ArrayXd& eta) : eta_(eta) {}
  LinearResponse operator+(const LinearResponse& rhs) const { return LinearResponse{eta_+rhs.get_eta()}; }
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, eta)
};
}
}

#endif

