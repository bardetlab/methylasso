#ifndef DESIGN_HPP
#define DESIGN_HPP

/*! \file Design.hpp
 * @brief implementation of Design class
 */

#include "typedefs.hpp"
#include "Binner.hpp"

namespace MethyLasso {
namespace params {

/// @brief Design represents a design matrix. Currently stored as a sparse matrix
class Design {
public:
  Design() =default;
  Design(const types::SpMat& X, const Eigen::ArrayXd& support) : X_(X), support_(support) {}
  
  unsigned get_nparams() const { return get_X().cols(); }
  
  METHLASSO_GET_CONST_OWN_DESCR(types::SpMat, X, dataset-level design matrix)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXd, support, support is a vector of size Nbins giving the value of the independent variable along these bins)
};

}
}

#endif
