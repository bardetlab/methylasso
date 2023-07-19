#include <RcppEigen.h>
using namespace Rcpp;

#include "signal_detection.hpp"
#include "methylation_helpers.hpp"
#include "difference_detection.hpp"
#include "weighted_1d_flsa.hpp"

RCPP_MODULE(MethyLasso_cpp) {
  using namespace Rcpp ;
  using namespace MethyLasso ;
  
  function("signal_detection_singlechr", &signal_detection_R,
           List::create(_["data"], _["optimize_lambda2"]=false, _["lambda2"]=50,
                        _["max_lambda2"]=1000, _["optimize_hyper"]=false, _["hyper"]=100, _["nouter"]=20, _["bf_per_kb"]=50,
                        _["tol_val"]=0.01, _["max_iter"]=30));
  
  function("signal_detection_fast_singlechr", &signal_detection_fast_R,
           List::create(_["data"], _["lambda2"]=50, _["tol_val"]=0.01, _["max_iter"]=1));
  
  function("segment_methylation_helper", &segment_methylation_helper, List::create(_["data"], _["tol_val"]=0.01));

  function("difference_detection_singlechr", &difference_detection_R,
           List::create(_["data"], _["ref"], _["optimize_lambda2"]=false,
                        _["lambda2"]=50, _["max_lambda2"]=1000, _["optimize_hyper"]=false, _["hyper"]=100, _["nouter"]=20,
                        _["bf_per_kb"]=50, _["tol_val"]=0.01, _["max_iter"]=30));
  
  function("weighted_1d_flsa", &weighted_1d_flsa, List::create(_["y"], _["wt"], _["lambda2"], _["lambda1"]=0));

  
}

