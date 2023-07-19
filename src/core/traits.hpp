#ifndef TRAITS_HPP
#define TRAITS_HPP

/*! \file traits.hpp
 * @brief Declares tags and traits used throughout
 */

#include <memory>

namespace MethyLasso {

//Method
struct Lasso {};
struct Mean {};

namespace binned {
class Binner; //forward declaration
}

namespace data {
class Observations; //forward declaration
class DataDesign; //
}

namespace params {
class Design; //forward declaration

/*! @brief Class which holds parameters which are optimized by a Fitter
 *
 * This class template must be specialized with a Method tag chosen for the
 * application. The specialization must have
 * - a public constructor, templated on a const ref Settings class as only
 *   argument, from which the initial values for all required parameters are to
 *   be extracted.
 * - any getters and setters for the aforementionned parameters. Use
 *   METHLASSO_GET_SET_OWN or related macros. There must be an ArrayXd parameter called
 *   beta (usually the main parameter vector) which will be used by binned::get_estimate.
 *   Use the _OWN macros to avoid useless copies of the data, Params is alive as
 *   long as Estimator.
 */
template<typename Method>
class Params;

/*! @brief Class which holds hyperparameters which are optimized by a Fitter
 *
 * This class template must be specialized with a Method tag chosen for the
 * application. The specialization must have
 * - a public constructor, templated on a const ref Settings class as only
 *   argument, from which the initial values for all required parameters are to
 *   be extracted.
 * - any getters and setters for the aforementionned hyperparameters. Use
 *   METHLASSO_GET_SET_DECL or related macros.
 * - a get_all() const function which returns all hyperparameters in an ArrayXd
 */
template<typename Method>
class Hyperparams;
}

namespace estimator {

/*! @brief Class to hold runtime-specific parameters for an Estimator
 */
template<typename Leg, typename Method>
struct Config;

/*! @brief Compile-time settings for Estimator
 */
template<typename Leg, typename Method>
struct EstimatorTraits;


/*! @brief function to create a (shared pointer to a) Binner from Observations, for a specific Leg and Method
 *
 * @param data the input data
 * @param conf the configuration options
 * @tparam Leg the leg tag of this design
 * @tparam Method the chosen method
 *
 * This function template must be specialized with an arbitrary Leg and Method tag
 * chosen for the application. It should create a Binner object, for example with the
 * help of binned::make_datasets_binner, and set remaining information accordingly by calling
 * - Binner::set_support
 */
template<typename Leg, typename Method>
std::shared_ptr<binned::Binner> make_binner_ptr(const data::Observations& data, const Config<Leg,Method>& conf);

/*! @brief function to create a Design from a Binner and a Config, for a specific Leg and Method
 *
 * @param binner the binner matrix, e.g. obtained through make_binner_ptr
 * @param conf the configuration options
 * @param data_design a single column of a data::DataDesign matrix
 * @tparam Leg the leg tag of this design
 * @tparam Method the chosen method
 *
 * This function template must be specialized with an arbitrary Leg and Method tag
 * chosen for the application. It should create a Design object, for example with the
 * help of params::make_signal_design or similar, using if necessary the information in Config
 */
template<typename Leg, typename Method>
params::Design make_design(const binned::Binner& binner, const Config<Leg,Method>& conf, const data::DataDesign& data_design);

/*! @brief Class which holds a specific implementation for a Fitter
 * 
 *  @tparam Method a tag for the specific method to use
 *  @tparam Library a method-specific library to use
 *
 * An instance is stored in the corresponding Estimator template.
 * The Fitter class template should be
 * specialized with a Method and a Library tag chosen for the application.
 * Typically, only partial specialization of Method is performed, leaving the
 * Library a template parameter. Estimator requires the Fitter to have the
 * following methods
 * - A public constructor which takes a Design class and a
 *   Config<Leg,Method> class, both by const ref. These are used to initialize
 *   any additional implementation-specific
 *   objects, such as the class to actually perform the fitting.
 * - An update_params(const params::Pseudodata& pseudo, bool damp = false, double new_param_weight = 1)
 *   call which performs the fit for these pseudodata
 * - A get_estimate() method which returns \f$\phi=X\beta\f$ as a binned::Estimate
 * - An update_hyperparams(bool damp = false, double new_param_weight = 1) call which 
 *   updates the hyperparameters if any, using damping if requested
 * - A get_hyperparams() method which returns them in an Eigen::ArrayXd
 * - A get_minus_log_prior() method which returns the -log(prior) as a double
 */
template<typename Method, typename Library>
class Fitter;
}
}


#endif

