#ifndef PSEUDODATA_HPP
#define PSEUDODATA_HPP

/*! \file Pseudodata.hpp
 * @brief Defines Pseudodata class
 */

#include <Eigen/Core>

#include "macros.hpp"

namespace MethyLasso {
namespace params {

//! @brief Pseudodata holds betahat and pseudoweight
/// @see \ref pseudodata_helpers.hpp
struct Pseudodata {
  Pseudodata() =default;
  Pseudodata(const Eigen::ArrayXd& betahat, const Eigen::ArrayXd& pseudoweight) :
     betahat_(betahat), pseudoweight_(pseudoweight) {}
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, betahat)
  METHLASSO_GET_SET_OWN(Eigen::ArrayXd, pseudoweight)
};
}
}

#endif

