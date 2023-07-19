#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

/*! \file typedefs.hpp
 * @brief simple typedefs
 */

#include <Eigen/Sparse>

namespace MethyLasso {
namespace types {

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, long> SpMat;

}
}

#endif
