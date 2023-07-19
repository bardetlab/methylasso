#ifndef OFFSET_HPP
#define OFFSET_HPP

#include <RcppEigen.h>
using namespace Rcpp;
#include <vector>
#include "typedefs.hpp"
#include "traits.hpp"
#include "macros.hpp"

#include "Estimator.hpp"
#include "fitter_mean.hpp"

namespace MethyLasso {

namespace params { class Design; } //forward declaration
namespace data { class Observations; }  //forward declaration
namespace binned { class Binner; }  //forward declaration

struct Offset {}; //declare tag

namespace estimator {

template<>
struct Config<Offset,Mean> {
  Config() {}
  Config(const Rcpp::List&) {}
};

//enable debug
template<>
struct EstimatorTraits<Offset,Mean> {
  //print debug info
  static const bool debug_summarizer = false;
  static const bool debug_fitter = false;
  static const bool debug_response = false;
  static const bool center_params = false;
  typedef void library;
};

template<>
std::shared_ptr<binned::Binner> make_binner_ptr<Offset,Mean>(const data::Observations& data, const Config<Offset,Mean>& conf);

template<>
params::Design make_design<Offset,Mean>(const binned::Binner& binner, const Config<Offset,Mean>& conf,
                                        const data::DataDesign& data_design);

} //namespace estimator

typedef estimator::Config<Offset,Mean> OffsetConfig;
typedef estimator::Estimator<Offset,Mean> OffsetEstimator;

}

#endif

