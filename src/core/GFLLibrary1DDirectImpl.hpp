#ifndef GFL_LIBRARY_1D_DIRECT_IMPL_HPP
#define GFL_LIBRARY_1D_DIRECT_IMPL_HPP

/*! \file GFLLibrary1DDirectImpl.hpp
 * @brief TODO
 */

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "GFLLibraryImpl.hpp"

namespace MethyLasso {
namespace base {

//policy class that implements fused lasso solution using the GFL library
//Direct implementation of Nick Johnson's dynamic programming algorithm 
//for exact O(n) calculation of the 1d fused lasso solution (at a given
//tuning parameter value).
class GFLLibrary1DDirectImpl : public GFLLibraryImpl {
public:
    
    //init by cold start
    GFLLibrary1DDirectImpl(int nrows) : GFLLibraryImpl(nrows, nrows), N_(nrows) {
        reset();
    }
    
    //cold start
    /*virtual*/ void reset();
    
    //return internal state
    /*virtual*/ GFLState_t get_state() const;
    
    //warm start
    /*virtual*/ void set_state(const GFLState_t& state);
    
    //run the optimization on the given data. The objective is
    // sum_i w_i(y_i-beta_i)^2 + lambda2 * sum_ij |beta_i-beta_j|
    // y, w and lambda2 are held constant, while beta starts at beta_init
    /*virtual*/ void optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge);
    
private:
    
    unsigned N_; //size of the fused lasso problem
    int ntrails_;
    std::vector<unsigned> trails_, breakpoints_;
    unsigned tsz_;
    
    std::vector<double> beta_impl_, x_, a_, b_, tm_, tp_;
    double *pbeta_, *px_, *pa_, *pb_, *ptm_, *ptp_; 
};

}
}

#endif

