#ifndef BINNED_RESIDUALS_HELPERS_HPP
#define BINNED_RESIDUALS_HELPERS_HPP

/*! \file binned_residuals_helpers.hpp
 * @brief functions to create BinnedResiduals objects
 */

#include "Binner.hpp"
#include "IRLSResiduals.hpp"
#include "BinnedResiduals.hpp"
#include "Estimate.hpp"

namespace MethyLasso {
namespace binned {

/// @brief Use Binner to build BinnedResiduals from IRLSResiduals
/*!
 *  @param z The IRLSResiduals object containing \f$z\f$ and \f$W\f$
 *  @param binner The binner \f$B_l\f$
 * 
 *  @return  a BinnedResiduals object, composed of zhat \f$\hat{z} \equiv \hat{W}_l^{-1} B_lWz\f$
 *  and binned weights \f$ \hat{W}_l \equiv B_l W B_l^\top \f$. Bins with no data will have a value of zero for
 *  the weight and zhat.
 *   
 */
BinnedResiduals bin_irls_residuals(const predictor::IRLSResiduals& z, const binned::Binner& binner);
}
}

#endif
