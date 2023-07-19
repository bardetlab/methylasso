#ifndef SIGNAL_DETECTION_HPP
#define SIGNAL_DETECTION_HPP

#include <Rcpp.h>

namespace MethyLasso {

Rcpp::List signal_detection_R(const Rcpp::DataFrame& data,
                              bool optimize_lambda2=false, double lambda2=50, double max_lambda2=1000,
                              bool optimize_hyper=false, double hyper=100,
                              unsigned nouter=20, double bf_per_kb=50, double tol_val=0.01,
                              unsigned max_iter=30);

Rcpp::List signal_detection_fast_R(const Rcpp::DataFrame& data, double lambda2=50,
                                   double tol_val=0.01, unsigned max_iter=1);



}

#endif

