#ifndef FUSED_LASSO_GAUSSIAN_ESTIMATOR_HPP
#define FUSED_LASSO_GAUSSIAN_ESTIMATOR_HPP

/*! \file FusedLassoGaussianEstimator.hpp
 * @brief TODO
 */

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "util.hpp"
#include "Settings.hpp"

namespace MethyLasso {
namespace base {

// A class that computes the 2D triangle grid fused lasso solution on some data.
// This class uses a gaussian model, hence assumes that weights are held constant.
// This class defines the full interface and implements the sparse part of the lasso,
// while the Library policy contains the dense implementation
template<typename Library>
class FusedLassoGaussianEstimator : public Library {
    
public:
    
    //initialize the problem with a triangle grid with nrows
    //requesting precision to be below a given convergence criterion
    //final beta value will be clamped if clamp>0
    FusedLassoGaussianEstimator(unsigned nrows, double tol, double max_lambda2)
       : Library(nrows), tol_(tol), max_lambda2_(max_lambda2),
         clamp_(Settings<FusedLassoGaussianEstimator<Library> >::get_clamp()) {}
    
    //approximate lasso solution computed on the first maxdiag counter-diagonals (included)
    //requesting precision to be below a given convergence criterion
    //final beta value will be clamped if clamp > 0
    //if maxdiag >= nrows, equivalent to the previous overload, computing the exact solution
    //for data beyond maxdiag, set beta to weighted average of all points
    FusedLassoGaussianEstimator(unsigned nrows, double tol, double max_lambda2, unsigned maxdiag)
       : Library(nrows, std::min(maxdiag,nrows)), tol_(tol), max_lambda2_(max_lambda2),
    clamp_(Settings<FusedLassoGaussianEstimator<Library> >::get_clamp()) {}
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2) {
        Library::optimize(y, w, lambda2, tol_);
    }

    //find the best lambda for the given betas
    double optimize_lambda(const std::vector<double>& beta) {
        double lambda = Library::optimize_lambda2(beta);
        lambda = std::min(max_lambda2_,lambda);
        return lambda;
    }
  
    /// @brief return the sum of minus log prior and hyperprior terms
    double get_minus_log_prior(const std::vector<double>& beta, double lambda) const {
      return Library::get_minus_log_prior(beta,lambda);
    }
  
    
    double get_tol() const { return tol_; }
    void set_tol(double tol) { tol_ = tol; }
    
    //return soft-thresholded value, corresponding to the problem
    // sum_i w_i(y_i-offset-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j| + lambda1 * sum_i |beta_i|
    //where the solution beta for a given lambda2 was already computed
    //values will be clamped if necessary
    //returns an empty vector if optimize has not been called
    std::vector<double> get(double offset=0., double lambda1=0.) const {
        return util::soft_threshold(clamp(Library::get_beta()), offset, lambda1);
    }
    
private:
    std::vector<double> clamp(std::vector<double> beta) const {
        //clamp values at +- clamp_ if needed
        if (clamp_>0) {
            for (double& i : beta) {
                i = std::min(clamp_, std::max(-clamp_, i));
            }
        }
        return beta;
    }
    
    double tol_;
    const double max_lambda2_;
    const double clamp_;
    
};

}
}

#endif

