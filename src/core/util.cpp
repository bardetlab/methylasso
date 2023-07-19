#include <vector>
#include <unsupported/Eigen/SpecialFunctions> // for log(Gamma)

#include "util.hpp"

namespace MethyLasso {
namespace util {

std::vector<double> soft_threshold(const std::vector<double>& beta,
                                   double eCprime, double lam1) {
    std::vector<double> phi;
    phi.reserve(beta.size());
    for (std::vector<double>::const_iterator it = beta.begin(); it != beta.end();
            ++it) {
        double val = *it - eCprime;
        phi.push_back((val > 0) ? std::max(0., val - lam1) : std::min(0., val + lam1));
    }
    return phi;
}

Eigen::ArrayXd log_beta(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y) {
  return x.lgamma() + y.lgamma() - (x+y).lgamma();
}



}
}



