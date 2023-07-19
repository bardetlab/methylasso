#ifndef BINNER_HELPERS_HPP
#define BINNER_HELPERS_HPP

/*! \file Binner.hpp
 * @brief Binner class and functions to create it it
 */

#include <memory>

#include "typedefs.hpp"
#include "macros.hpp"
#include "Binner.hpp"

namespace MethyLasso {
namespace binned {

/// @brief Create BinDesign for evenly sized bins with labels from 0 to N
/*! Contains final boundary (size is N+1) */
std::shared_ptr<BinDesign> make_even_bins(double lower, double upper, unsigned Nbins);

/// @brief places each of the integer-valued elements in data in a separate bin, in increasing order. Two equal values are in the same bin
std::shared_ptr<BinDesign> make_unique_bins(const Eigen::ArrayXi& pos);

/// @brief bin data given a boundary map in the form [a,b) [b,c) ... [y,z] with labels from 0 to N-1
/*! Return in Nbins x Ndata sparse matrix form.
 * by default, drops unused bins (rows with only zeroes)
 */
std::shared_ptr<Binner> make_binner_using_bin_design(const Eigen::ArrayXd& data, const std::shared_ptr<BinDesign>& design, bool drop = true);

/// @brief split the array in data into groups defined by dset, which must start at 1 without gaps
std::vector<Eigen::ArrayXd> split_by_dataset(const Eigen::ArrayXi& dset, const Eigen::ArrayXd& data);

/// @brief call a funnction on each dataset indexed by the dset vector, and build block diagonal binner on the result
/*!
 * @param pos_by_dset a vector of vectors of positions for each dataset, as returned by split_by_dataset
 * @param make_binner a function that takes a vector of positions as input, and produces a std::shared_ptr<Binner>
 * @param args any additional arguments to make_binner, if any
 * 
 * @note could be implemented using std::accumulate, but its definition forbids modifying any of
 *  the arguments to binary_op. Here, we need to, to avoid making copies of the Binner object
 */
template<class BinnerMaker, class... BinnerMakerArgs>
std::shared_ptr<binned::Binner> make_datasets_binner(const std::vector<Eigen::ArrayXd>& pos_by_dset, BinnerMaker&& make_binner, BinnerMakerArgs&&... args) {
  //bin each dataset and combine with previous one
  std::shared_ptr<binned::Binner> ret(make_binner(pos_by_dset[0], std::forward<BinnerMakerArgs>(args)...));
  for (unsigned i=1; i<pos_by_dset.size(); ++i) {
    auto new_b = make_binner(pos_by_dset[i], std::forward<BinnerMakerArgs>(args)...);
    ret->extend_diagonal(*new_b);
  }
  return ret;
}


}
}

#endif
