#ifndef WEIGHTED_1D_FLSA_HPP
#define WEIGHTED_1D_FLSA_HPP

#include <RcppEigen.h>

namespace MethyLasso {

Eigen::VectorXd weighted_1d_flsa(const Eigen::VectorXd& y, const Eigen::VectorXd& wt, double lambda2, double lambda1 = 0);

}

#endif

