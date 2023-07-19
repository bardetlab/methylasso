#ifndef MEAN_VECTOR_HPP
#define MEAN_VECTOR_HPP

/*! \file MeanVector.hpp
 * @brief Defines MeanVector class
 */

#include <Eigen/Core>

#include "macros.hpp"

namespace MethyLasso {
namespace predictor {

//! @brief MeanVector holds \f$\mu = g(\eta)\f$
struct MeanVector {
  MeanVector() =default;
  MeanVector(const Eigen::ArrayXd& mu) : mu_(mu) {}
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, mu)
};
}
}

#endif

