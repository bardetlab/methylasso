#ifndef DATA_DESIGN_HELPERS_HPP
#define DATA_DESIGN_HELPERS_HPP

/*! \file data_design_helpers.hpp
 * @brief functions to make dataset-level design matrices
 */

#include <Eigen/Core>

#include "DataDesign.hpp"
#include "Observations.hpp"

namespace MethyLasso {
namespace data {

/*! @brief build design matrix for signal calculation
 */
DataDesign build_signal_design(const Observations& obs);

/*! @brief build design matrix for difference calculation
 */
DataDesign build_difference_design(const Observations& obs, unsigned ref_condition = 1);

}
}

#endif
