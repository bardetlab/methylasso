#ifndef GFL_LIBRARY_1D_ADAMTS_IMPL_HPP
#define GFL_LIBRARY_1D_ADAMTS_IMPL_HPP

/*! \file GFLLibrary1DAdamtsImpl.hpp
 * @brief TODO
 */

#include <Rcpp.h>
using namespace Rcpp;
#include <vector>

#include "GFLLibraryImpl.hpp"

namespace MethyLasso {
namespace base {

std::vector<std::vector<int> > linear_chain(int nrows);

//policy class that implements fused lasso solution using the GFL library
class GFLLibrary1DAdamtsImpl : public GFLLibraryImpl {
public:
    
    //init by cold start
    GFLLibrary1DAdamtsImpl(int nrows) : GFLLibraryImpl(nrows, nrows), N_(nrows) {
        store_trails(nrows);
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
    
    void store_trails(int nrows);
    
    unsigned N_; //size of the fused lasso problem
    int ntrails_;
    std::vector<unsigned> trails_, breakpoints_;
    unsigned tsz_;
    
    std::vector<double> beta_impl_, z_, u_;
};

}
}

#endif

