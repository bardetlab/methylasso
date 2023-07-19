#include <Rcpp.h>
#include <vector>
#include <numeric> //std::iota
#include <RcppGSL.h>

#include "GFLLibrary1DAdamtsImpl.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm


namespace MethyLasso {
namespace base {


std::vector<std::vector<int> > linear_chain(int nrows) {
    std::vector<std::vector<int> > chains;
    std::vector<int> current(nrows);
    std::iota(std::begin(current),std::end(current),0); //current : 0, 1, ..., nrows-1
    chains.push_back(current);
    return(chains);
}

void GFLLibrary1DAdamtsImpl::store_trails(int nrows) {
    const std::vector<std::vector<int> > chains = linear_chain(nrows);
    //
    trails_.clear();
    breakpoints_.clear();
    for (std::vector<std::vector<int> >::const_iterator it = chains.begin() ;
         it != chains.end(); ++it) {
        if (trails_.size()>0) breakpoints_.push_back(trails_.size());
        trails_.insert(trails_.end(), it->begin(), it->end());
    }
    if (trails_.size()>0) breakpoints_.push_back(trails_.size());
    ntrails_ = breakpoints_.size();
    tsz_ = trails_.size();
}

void GFLLibrary1DAdamtsImpl::reset() {
    //setup initial values for a cold start
    set_ninner(0);
    beta_impl_ = std::vector<double>(N_,0);
    z_ = std::vector<double>(tsz_,0);
    u_ = std::vector<double>(tsz_,0);
}

void GFLLibrary1DAdamtsImpl::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
    //perform optimization on the C side
    double* py = const_cast<double*>(&y[0]);
    double* pw = const_cast<double*>(&w[0]);
    double alpha = get_alpha();
    int counter = graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha, get_inflate(), get_ninner_max(), converge,
                                               &beta_impl_[0], &z_[0], &u_[0]);
    set_alpha(alpha);
    set_beta(beta_impl_);
    set_ninner(get_ninner() + counter);
}

GFLLibrary1DAdamtsImpl::GFLState_t GFLLibrary1DAdamtsImpl::get_state() const {
    return Rcpp::List::create(_["z"]=z_, _["u"]=u_,  _["alpha"]=get_alpha(),  _["beta"]=beta_impl_,
                              _["counter"]=get_ninner());
}

void GFLLibrary1DAdamtsImpl::set_state(const GFLState_t& state) {
    if (state.containsElementNamed("u") && Rcpp::as<std::vector<double> >(state["u"]).size() == tsz_
          && state.containsElementNamed("beta") && state.containsElementNamed("alpha")
          && state.containsElementNamed("counter")) {
        z_ = Rcpp::as<std::vector<double> >(state["z"]);
        u_ = Rcpp::as<std::vector<double> >(state["u"]);
        beta_impl_ = Rcpp::as<std::vector<double> >(state["beta"]);
        set_alpha(Rcpp::as<double>(state["alpha"]));
        set_ninner(Rcpp::as<unsigned>(state["counter"]));
    }
}

}
}

