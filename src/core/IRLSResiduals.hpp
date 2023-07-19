#ifndef IRLS_RESIDUALS_HPP
#define IRLS_RESIDUALS_HPP

/*! \file IRLSResiduals.hpp
 */

#include <Eigen/Core>
#include "macros.hpp"

namespace MethyLasso {
namespace predictor {

struct IRLSResiduals {
public:
    IRLSResiduals(const Eigen::ArrayXd& residuals, const Eigen::ArrayXd& weights)
        : residuals_(residuals), weights_(weights) {}
    METHLASSO_GET_SET_OWN(Eigen::ArrayXd, residuals)
    METHLASSO_GET_SET_OWN(Eigen::ArrayXd, weights)
};

}
}

#endif

