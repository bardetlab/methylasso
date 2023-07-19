#ifndef GFL_LIBRARY_IMPL_HPP
#define GFL_LIBRARY_IMPL_HPP

/*! \file GFLLibraryImpl.hpp
 * @brief TODO
 */

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "Settings.hpp"

namespace MethyLasso {
namespace base {

//policy class that implements fused lasso solution using the GFL library
class GFLLibraryImpl {
  public:
    
    GFLLibraryImpl(unsigned nrows, unsigned maxdiag) : counter_(0),
    inflate_(Settings<GFLLibraryImpl>::get_inflate()), ninner_max_(Settings<GFLLibraryImpl>::get_ninner_max()),
    alpha_(Settings<GFLLibraryImpl>::get_alpha()) {}
 
  typedef Rcpp::List GFLState_t;
    
  //cold start, implementation-specific
  virtual void reset() =0;
  
  //return internal state
  virtual GFLState_t get_state() const =0;
  
  //warm start
  virtual void set_state(const GFLState_t& state) =0;
  
  //run the optimization on the given data. The objective is
  // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
  // y, w and lambda2 are held constant, while beta starts at beta_init
  virtual void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) =0;
  
  //return modified quantities
  std::vector<double> get_beta() const {
    return beta_;
  }
  
  unsigned get_ninner() const {
    return counter_;
  }
  
  unsigned get_ninner_max() const {
    return ninner_max_;
  }
  
  double get_alpha() const {
    return alpha_;
  }
 
  double get_inflate() const {
    return inflate_;
  }
  
protected:
  void set_beta(const std::vector<double>& beta) { beta_ = beta; }
  
  void set_ninner(unsigned ninner) { counter_ = ninner; }
  
  void set_ninner_max(unsigned ninner_max) { ninner_max_ = ninner_max; }
  
  void set_alpha(double alpha) { alpha_ = alpha; }
  
  void set_inflate(double inflate) { inflate_ = inflate; }
  
  private:
    unsigned counter_;
    double inflate_;
    int ninner_max_;
    double alpha_;
    std::vector<double> beta_;
};

}
}

#endif

