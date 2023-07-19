#include <Rcpp.h>
#include <vector>
#include <numeric> //std::iota
#include <RcppGSL.h>

#include "GFLLibrary1DDirectImpl.hpp"
#include "gfl_tf.h" // tf_dp_weight

namespace MethyLasso {
namespace base {

void GFLLibrary1DDirectImpl::reset() {
    //setup initial values for a cold start
    set_ninner(0);
    
    beta_impl_ = std::vector<double>(N_,0);
    pbeta_ = &beta_impl_[0];
    
    x_ = std::vector<double>(2*N_,0);
    a_ = std::vector<double>(2*N_,0);
    b_ = std::vector<double>(2*N_,0);
    px_ = &x_[0];
    pa_ = &a_[0];
    pb_ = &b_[0];
    
    tm_ = std::vector<double>(N_-1,0);
    tp_ = std::vector<double>(N_-1,0);
    ptm_ = &tm_[0];
    ptp_ = &tp_[0];
}

void GFLLibrary1DDirectImpl::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
  
  //some safety checks
  if (!(N_ == y.size() && N_ == w.size())) throw std::invalid_argument("Invalid sizes for GFLLibrary1DDirectImpl inputs! Expected " + std::to_string(N_) + " got " + std::to_string(y.size()));
  
  //perform optimization on the C side
  double* py = const_cast<double*>(&y[0]);
  double* pw = const_cast<double*>(&w[0]);
  
  tf_dp_weight(N_, py, pw, lambda2, pbeta_, px_, pa_, pb_, ptm_, ptp_);
  
  set_beta(beta_impl_);
  //Rcpp::Rcout << "y[33]= " << y[33] << " w[33]= " << w[33] << " beta[33]= " << beta_impl_[33] << " \n";
  set_ninner(get_ninner() + 1);
}

GFLLibrary1DDirectImpl::GFLState_t GFLLibrary1DDirectImpl::get_state() const {
    return Rcpp::List::create();
}

void GFLLibrary1DDirectImpl::set_state(const GFLState_t&) {}

}
}


