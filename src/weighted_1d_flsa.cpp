#include "weighted_1d_flsa.hpp"
#include "FusedLassoGaussianEstimator.hpp"
#include "GFLLibrary1D.hpp"

namespace MethyLasso {

Eigen::VectorXd weighted_1d_flsa(const Eigen::VectorXd& y, const Eigen::VectorXd& wt, double lambda2, double lambda1) {
    const unsigned ndata = y.size();
    if (wt.size() != ndata) Rcpp::stop("Input vectors must be of same size!");
    const double tol_val = 1e-2; //unused in 1D case
    const double max_lambda2 = lambda2+1; //unused in fixed lambda case
    base::FusedLassoGaussianEstimator<base::GFLLibrary1D> flo(ndata, tol_val, max_lambda2);

    std::vector<double> y_std(y.data(), y.data() + ndata);
    std::vector<double> wt_std(wt.data(), wt.data() + ndata);
    flo.optimize(y_std,wt_std,lambda2);
    std::vector<double> beta_std(flo.get(0, lambda1));
    Eigen::Map<Eigen::ArrayXd> beta(beta_std.data(), beta_std.size());
    
    return Eigen::VectorXd(beta);
}

}

