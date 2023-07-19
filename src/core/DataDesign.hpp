#ifndef DATA_DESIGN_HPP
#define DATA_DESIGN_HPP

/*! \file DataDesign.hpp
 * @brief implementation of DataDesign class and friend functions
 */

#include "typedefs.hpp"
#include "macros.hpp"

namespace MethyLasso {
namespace data {

/// @brief DataDesign represents a design matrix at dataset-level. Currently stored as a sparse matrix
class DataDesign {
public:
  DataDesign() =default;
  DataDesign(const types::SpMat& X) : X_(X) {}
  
  unsigned get_num_estimators() const { return get_X().cols(); }
  Eigen::SparseMatrix<double> get_col(unsigned i) const { return get_X().col(i); }
  
  METHLASSO_GET_CONST_OWN_DESCR(types::SpMat, X, dataset-level design matrix)
};

}
}

#endif
