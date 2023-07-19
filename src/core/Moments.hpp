#ifndef MOMENTS_HPP
#define MOMENTS_HPP

/*! \file Moments.hpp
 * @brief Defines Moments class
 */

#include <Eigen/Core>

#include "macros.hpp"
#include "MeanVector.hpp"

namespace MethyLasso {
namespace predictor {

//! @brief Moments holds the first two moments of a distribution (\f$\mu\f$ and \f$V\f$)
/*! @note the variance is always with respect to the parameter \f$\mu\f$ to be estimated, and
 *    might be different from the variance of the distribution itself. This happens in particular if
 *    the weight of a given observation is different from 1, which affects the variance of the
 *    estimated \f$\mu\f$.
 */
struct Moments {
  Moments(const Eigen::ArrayXd& mean, const Eigen::ArrayXd& variance) : mean_(mean), variance_(variance) {}
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, mean)
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, variance)
};
}
}

#endif

