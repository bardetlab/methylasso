#ifndef ESTIMATE_HELPERS_HPP
#define ESTIMATE_HELPERS_HPP

/*! \file estimate_helpers.hpp
 * @brief functions to create Estimate objects
 */

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Design.hpp"
#include "Estimate.hpp"

namespace MethyLasso {
namespace binned {

class Binner;
class BinnedResiduals;

/// @brief Use Design to build Estimate from Parameters
/*!
 *  @param params The Parameters object containing \f$\beta\f$
 *  @param X The design matrix \f$X\f$
 * 
 *  @return \f$\phi \equiv X\beta\f$
 *   
 */
template<typename Parameters>
Estimate get_estimate(const Parameters& params, const params::Design& des) {
  auto X = des.get_X();
  return Estimate{(X*params.get_beta().matrix()).array()};
}

}
}

#endif
