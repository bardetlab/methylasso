#ifndef PREDICTOR_HELPERS_HPP
#define PREDICTOR_HELPERS_HPP

/*! \file predictor_helpers.hpp
 * @brief functions to create LinearResponse and MeanVector objects
 */

#include "Binner.hpp"
#include "Estimate.hpp"
#include "LinearResponse.hpp"
#include "MeanVector.hpp"
#include "Moments.hpp"
#include "IRLSResiduals.hpp"
#include "Observations.hpp"
#include "link_functions.hpp"

namespace MethyLasso {
namespace predictor {

/// @brief Build LinearResponse from a single Estimate and Binner pair
/*!
 *  @param est The Estimate object containing \f$\phi\f$
 *  @param B The binner matrix \f$X\f$
 * 
 *  @return \f$\eta \equiv B^\top \phi\f$
 *  @note To get \f$\sum_i \eta_i\f$ for a series of estimates and associated binners,
 *  simply sum the LinearResponse objects.
 *   
 */
LinearResponse get_linear_response(const binned::Estimate& est, const binned::Binner& B);

/// @brief Build Moments from a MeanVector and a distribution
/*!
 *  @param mean the MeanVector \f$\mu\f$
 *  @param Nobs the number of observations (or weight) for that mean
 *  @param dist the Distribution object considered. Will use hyperparameter value if relevant.
 *  @tparam Distribution The type of distribution, deduced automatically
 * 
 *  @return \f$\mu \equiv g^{-1}(\eta)\f$
 *   
 */
template<class Distribution>
Moments get_moments_from_distribution(const MeanVector& mean, const data::Observations& obs, const Distribution& dist) {
  return Moments{mean.get_mu(),dist.get_variance(obs.get_Nobs(),mean.get_mu())};
}

/// @brief Return maximum absolute deviation between the two mean vectors
double get_max_absolute_deviation(const MeanVector& mu1, const MeanVector& mu2);


/// @brief Return minus log pdf of obs given mu and distribution dist (which can contain hyperparameters)
template<typename Distribution>
double get_minus_log_pdf(const data::Observations& obs, const MeanVector& mu, const Distribution& dist) {
  return dist.get_minus_log_pdf(obs.get_observed(), obs.get_Nobs(), mu.get_mu());
}

}
}

#endif
