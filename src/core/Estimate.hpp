#ifndef ESTIMATE_HPP
#define ESTIMATE_HPP

/*! \file Estimate.hpp
 * @brief Defines Estimate class
 */

#include <Eigen/Core>

#include "macros.hpp"

namespace MethyLasso {
namespace binned {

//! @brief Estimate holds \f$\phi\f$
struct Estimate {
  Estimate(const Eigen::ArrayXd& phi) : phi_(phi) {}
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, phi)
};
}
}

#endif

