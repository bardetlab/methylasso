#ifndef ESTIMATOR_TRAITS_HPP
#define ESTIMATOR_TRAITS_HPP

/*! \file EstimatorTraits.hpp
 * @brief Defines the EstimatorTraits class
 */

namespace MethyLasso {
namespace estimator {
//! @brief EstimatorTraits holds compile-time options
/*!
 * If an implementation of the estimator needs to change this variable,
 * template specialization allows to do it case-by-case.
 * 
 * @tparam Leg the Leg of this estimator
 * @tparam Method the prior placed on the parameters to be estimated
 */
template<typename Leg, typename Method>
struct EstimatorTraits {
  //this flag should be always false, but can be overridden in a specific template specialization if needed
  static const bool debug_summarizer = false;
  static const bool debug_fitter = false;
  static const bool debug_response = false;
  //which class to use to do the actual fitting
  typedef void library;
};
}
}

#endif

