#ifndef UNIQUE_POSITIONS_HPP
#define UNIQUE_POSITIONS_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "typedefs.hpp"
#include "traits.hpp"
#include "macros.hpp"

#include "Estimator.hpp"
#include "EstimatorGroup.hpp"
#include "GFLLibrary1D.hpp"
#include "fitter_lasso.hpp"

namespace MethyLasso {

namespace params { class Design; } //forward declaration
namespace data { class Observations; }  //forward declaration
namespace binned { class Binner; }  //forward declaration

struct UniquePositions {}; //declare tag

namespace estimator {

template<>
struct Config<UniquePositions,Lasso> {
  Config(double lambda2=1, double max_lambda2=1000, double lambda1=0, unsigned max_iter=20,
         double tol_val=0.01, unsigned nouter=20) :
     lambda2_(lambda2), max_lambda2_(max_lambda2), lambda1_(lambda1), max_iter_(max_iter),
     tol_val_(tol_val), nouter_(nouter) {}
  
  Config(const Rcpp::List& conf) :
     lambda2_(conf["lambda2"]), max_lambda2_(conf["max_lambda2"]), lambda1_(conf["lambda1"]),
     max_iter_(conf["max_iter"]), tol_val_(conf["tol_val"]), nouter_(conf["nouter"]) {}
  
  METHLASSO_GET_CONST_SIMPLE(double, lambda2)
  METHLASSO_GET_CONST_SIMPLE(double, max_lambda2)
  METHLASSO_GET_CONST_SIMPLE(double, lambda1)
  METHLASSO_GET_CONST_SIMPLE(unsigned, max_iter)
  METHLASSO_GET_CONST_SIMPLE(double, tol_val)
  METHLASSO_GET_CONST_SIMPLE(unsigned, nouter)
};

//enable debug
template<>
struct EstimatorTraits<UniquePositions,Lasso> {
  //print debug info
  static const bool debug_summarizer = false;
  static const bool debug_fitter = false;
  static const bool debug_response = false;
  static const bool center_params = true;
  typedef base::GFLLibrary1D library;
};

/// in UniquePositions, make_binner_ptr will build a binner which maps each position to a bin, in order
template<>
std::shared_ptr<binned::Binner> make_binner_ptr<UniquePositions,Lasso>(const data::Observations& data,
                                                                    const Config<UniquePositions,Lasso>& conf);

/// in UniquePositions, make_design will build a design based on the groups, assuming each dataset is binned separately
template<>
params::Design make_design<UniquePositions,Lasso>(const binned::Binner& binner, const Config<UniquePositions,Lasso>& conf,
                                               const data::DataDesign& data_design);

} //namespace estimator

typedef estimator::Config<UniquePositions,Lasso> UniquePositionsConfig;
typedef estimator::Estimator<UniquePositions,Lasso> UniquePositionsEstimator;
typedef estimator::EstimatorGroup<UniquePositions,Lasso> UniquePositionsEstimatorGroup;
          

}

#endif

