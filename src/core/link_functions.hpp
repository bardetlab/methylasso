#ifndef LINK_FUNCTIONS_HPP
#define LINK_FUNCTIONS_HPP

/*! \file link_functions.hpp
 * @brief Implements several link functions
 *
 * Currently, IdentityLink, LogLink and LogitLink are implemented.
 */

#include <Eigen/Core>

#include "IRLSResiduals.hpp"
#include "macros.hpp"
#include "LinearResponse.hpp"
#include "MeanVector.hpp"
#include "Moments.hpp"
#include "Observations.hpp"
#include "traits.hpp"

namespace MethyLasso {
namespace link {

/*! @defgroup link_functions All link functions specialize the following templates
 *  @{
 */

/// @brief evaluates the link function \f$\eta = g(\mu)\f$
template<typename LinkName>
predictor::LinearResponse get_eta(const predictor::MeanVector&);

/// @brief evaluates the inverse of the link function \f$\mu = g^{-1}(\eta)\f$
template<typename LinkName>
predictor::MeanVector get_mu(const predictor::LinearResponse&);

/// @brief Restrict values of input vector to the support of g.
/*! For example, in the logit, this would be \f$]-1,1[\f$. This is only called for the initial guess. */
template<typename LinkName>
predictor::MeanVector restrict_to_support(const predictor::MeanVector&);

/// @brief Given an initial value for \f$\mu^0\f$ and \f$V^0\f$, return \f$\hat{\phi}^0\f$ and \f$W^0\f$.
/*!
 * Usually, predictor::IRLSResiduals contains \f$(z,W)\f$, but here
 *    we return \f$(\hat{\phi},W)\f$, which will subsequently be passed to the estimator, which assumes \f$\beta=0\f$
 *    for the initialization
*/
template<typename LinkName>
predictor::IRLSResiduals transform_initial_guess(const predictor::Moments&);

/// @brief builds the next set of IRLS residuals from the provided moments and data
/*!
 *  @param obs the actual observations
 *  @param moments The Moments object containing \f$\mu\f$ and \f$V\f$
 * 
 *  @return IRLSResiduals object, containing \f$z \equiv g'(\mu)(d-\mu)\f$
 *  and \f$W \equiv \frac{1}{V g'^2(\mu)}\f$
 */
template<typename LinkName>
predictor::IRLSResiduals get_irls_residuals(const data::Observations&, const predictor::Moments&);

/// @}

/*! @defgroup link_functions_impl Link functions currently implemented
 *  @{
 */
/// @brief \f$g(x) = x\f$
class Identity {};

/// @brief \f$g(x) = \log x \quad g^{-1}(x) = \exp x\f$
class Logarithmic {};

/// @brief \f$g(x) = \log \frac{x}{1-x} \quad g^{-1}(x) = \frac{1}{e^{-x}+1}\f$
class Logit {};

/// @}

}
}

#endif

