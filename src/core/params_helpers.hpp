#ifndef PARAMS_HELPERS_HPP
#define PARAMS_HELPERS_HPP

/*! \file params_helpers.hpp
 * @brief functions to manipulate Params objects
 */

#include "traits.hpp"
#include "Pseudodata.hpp"

namespace MethyLasso {
namespace params {

/// @brief Center a Params vector by subtracting its weighted average
/*!
 *  @param par the Params to be centered
 *  @param pseudo the Pseudodata which contain the weights to be used
 * 
 *  @return \f$\beta\f$ such that \f$1^\top \tilde{W}_l^k \beta = 0\f$
 *   
 */
template<typename Method>
Params<Method> center_params(const Params<Method>& par, const Pseudodata& pseudo) {
  //extract data
  const Eigen::ArrayXd weight = pseudo.get_pseudoweight();
  const Eigen::ArrayXd beta = par.get_beta();
  //compute and subtract mean
  double mean = (beta*weight).sum()/weight.sum();
  const Eigen::ArrayXd new_beta = beta - mean;
  return Params<Method>{new_beta};
}


}
}

#endif
