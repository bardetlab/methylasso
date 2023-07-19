#ifndef PSEUDODATA_HELPERS_HPP
#define PSEUDODATA_HELPERS_HPP

/*! \file pseudodata_helpers.hpp
 * @brief functions to create Pseudodata objects
 */

#include "Design.hpp"
#include "BinnedResiduals.hpp"
#include "Pseudodata.hpp"

namespace MethyLasso {
namespace params {

/// @brief Use Design to build Pseudodata from BinnedResiduals
/*!
 *  @param s The BinnedResiduals object containing \f$\hat{\phi}\f$ and \f$W_l\f$
 *  @param des The design matrix \f$X\f$
 *  @param par the previous value of the parameters
 * 
 *  @return  a Pseudodata object, composed of betahat \f$\hat{\beta} \equiv \tilde{W}^{-1} X^\top W_l \hat{z} + \beta\f$
 *  and pseudoweights \f$ \tilde{W} \equiv X^\top W_l X \f$. Bins with no data will have a value of zero for
 *  the pseudoweight and beta for betahat (amounts to using the pseudoinverse for \f$\tilde{W}\f$).
 *   
 */
template<typename Params>
Pseudodata design_pseudodata(const binned::BinnedResiduals& s, const params::Design& des, const Params& par) {
  auto X = des.get_X();
  //compute pseudoweight
  Eigen::SparseMatrix<double> xwx = X.transpose()*s.get_weight().matrix().asDiagonal()*X;
  Eigen::ArrayXd pseudoweight = xwx.diagonal().array();
  //compute betahat
  Eigen::ArrayXd wt_inv = pseudoweight.inverse();
  wt_inv = wt_inv.isFinite().select(wt_inv,Eigen::ArrayXd::Zero(wt_inv.size())); //use pseudoinverse
  Eigen::ArrayXd betahat = ( (wt_inv.matrix().asDiagonal() * X.transpose())
                             * (s.get_weight().matrix().asDiagonal() * s.get_zhat().matrix()) ).array()
                           + par.get_beta();
  return Pseudodata{betahat,pseudoweight};
}

}
}

#endif
