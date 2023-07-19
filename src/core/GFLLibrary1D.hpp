#ifndef GFL_LIBRARY_1D_HPP
#define GFL_LIBRARY_1D_HPP

/*! \file GFLLibrary1D.hpp
 * @brief TODO
 */

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>
#include <memory>
#include <cmath>

#include "Settings.hpp"
#include "GFLLibraryImpl.hpp"
#include "GFLLibrary1DDirectImpl.hpp"

namespace MethyLasso {
namespace base {

//policy class that implements 1D fused lasso solution using the GFL library
class GFLLibrary1D {
public:
    
    typedef Rcpp::List GFLState_t;
    
    //init by cold start
    GFLLibrary1D(unsigned nrows) {
        strategy_ = std::make_shared<GFLLibrary1DDirectImpl>(nrows);
    }

    //reset
    void reset() { strategy_->reset(); }
    
    //return internal state
    GFLState_t get_state() const { return strategy_->get_state(); }
    
    //warm start
    void set_state(const GFLState_t& state) { strategy_->set_state(state); }
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
      strategy_->optimize(y,w,lambda2,converge);
    }
    
    //optimize lambda2
    //1D setting has analytical formula
    //lambda = (N-2+alphatilde) /( sum_i |beta_i-beta_i+1| + betatilde)
    //here alphatilde=0 and betatilde=0 (Jeffreys prior)
    double optimize_lambda2(const std::vector<double>& beta) {
      double lambda=0;
      for (unsigned i=1; i<beta.size(); ++i) lambda += std::abs(beta[i]-beta[i-1]);
      const double alphatilde=0, betatilde=0;
      lambda = (beta.size()-2+alphatilde)/(lambda+betatilde);
      return lambda;
    }
    
    /// @brief return the sum of minus log prior and hyperprior terms
    /*! We are using 1D fused lasso with Jeffreys prior (it is invariant by global translation) */
    double get_minus_log_prior(const std::vector<double>& beta, double lambda) const {
      double ret=0;
      for (unsigned i=1; i<beta.size(); ++i) ret += std::abs(beta[i]-beta[i-1]);
      ret = lambda*ret - (beta.size()-2)*std::log(lambda) + (beta.size()-1)*std::log(2);
      return ret;
    }
    
    //return modified quantities
    std::vector<double> get_beta() const { return strategy_->get_beta(); }
    
    unsigned get_ninner() const { return strategy_->get_ninner(); }

    double get_alpha() const { return strategy_->get_alpha(); }
    
//protected:
    //to avoid direct destruction by user (follows policy-based class design)
    ~GFLLibrary1D() = default;
    
private:
    std::shared_ptr<GFLLibraryImpl> strategy_;
};

}
}

#endif

