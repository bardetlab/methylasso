#ifndef UTIL_HPP
#define UTIL_HPP

/*! \file util.hpp
 * @brief Miscellaneous utility functions
 *
 * Currently, has soft_threshold, and log_beta
 */

#include <Rcpp.h>
#include <vector>

namespace MethyLasso {
namespace util {

/// @brief apply soft-thresholding operator
std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1);

/// @brief log(beta(x,y)) function
Eigen::ArrayXd log_beta(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y);

}
}

#endif

