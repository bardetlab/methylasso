#ifndef DESIGN_HELPERS_HPP
#define DESIGN_HELPERS_HPP

/*! \file design_helpers.hpp
 * @brief implementation of functions to create Design objects
 */

#include "typedefs.hpp"

namespace MethyLasso {
namespace data {
class DataDesign;
}
namespace binned {
class Binner;
}
namespace params {
class Design;

/// @brief Simple identity design
///
/// If binner has N bins (including empty bins), Design will be N x N
Design make_identity_design(const binned::Binner& binner);

/// @brief design by groups
/*! 
 * @param the Binner matrix class used to map data to bins (contains dataset info)
 * @param grouped_design the way each dataset group in the Binner is related to each other
 * 
 * Suppose the Binner matrix has N rows, and a total of K bins per dataset (counting unused ones).
 * The design matrix will be N x K. The coefficient of a row belonging to dataset group i will be grouped_design(i)
 */
Design make_group_design(const binned::Binner& binner, const data::DataDesign& grouped_design);

}
}

#endif
