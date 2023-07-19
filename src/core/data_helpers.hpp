#ifndef DATA_HELPERS_HPP
#define DATA_HELPERS_HPP

/*! \file data_helpers.hpp
 * @brief functions to group datasets
 */

#include <Eigen/Core>

#include "DataDesign.hpp"

namespace MethyLasso {
namespace data {

/*! @brief group datasets based on a specified design
 * 
 * All datasets which have the same coefficient value get the same group ID. Zero coefficients are discarded.
 * 
 * @param data_design a single column of a data::DataDesign matrix
 * @param tol data_design coefficients which are closer than tol are considered equal,
 * and corresponding datasets are therefore grouped.
 */
Eigen::ArrayXi group_datasets(const DataDesign& data_design, double tol=1e-6);

/*! @brief transform a data_design object into the corresponding design for a group of datasets
 * 
 * @param data_design a single column of a data::DataDesign matrix
 * @param tol data_design coefficients which are closer than tol are considered equal,
 * and corresponding datasets are therefore grouped.
 */
DataDesign get_group_design(const DataDesign& data_design, const Eigen::ArrayXi& group_ids);


}
}

#endif
