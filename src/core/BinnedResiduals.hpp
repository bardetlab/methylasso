#ifndef BINNED_RESIDUALS_HPP
#define BINNED_RESIDUALS_HPP

/*! \file BinnedResiduals.hpp
 * @brief Defines BinnedResiduals class
 */

#include <Eigen/Core>

#include "macros.hpp"

namespace MethyLasso {
namespace binned {

//! @brief BinnedResiduals holds zhat and weight
class BinnedResiduals {
    public:
  BinnedResiduals() =default;
  BinnedResiduals(const Eigen::ArrayXd& zhat, const Eigen::ArrayXd& weight) :
     zhat_(zhat), weight_(weight) {}
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, zhat)
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, weight)
};
}
}

#endif

