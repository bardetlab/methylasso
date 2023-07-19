#ifndef BINNER_HPP
#define BINNER_HPP

#include <Rcpp.h>
#include <memory>


/*! \file Binner.hpp
 * @brief implementation of Binner class
 */

#include "typedefs.hpp"
#include "macros.hpp"
#include "BinDesign.hpp"

namespace MethyLasso {
namespace binned {
/// @brief Binner class represents a binner matrix
/*! A binner matrix is a nbins x ndta matrix \f$B\f$ such that \f$b_{ij}=1\f$ if data point \f$j\f$ belongs to bin \f$i\f$.
 * Currently, we use std::map to store the Binner, and not a sparse matrix.
 */
class Binner {
public:

  /// @brief build the Binner from a sparse matrix representing a single dataset
  Binner(const types::SpMat& binner, const std::shared_ptr<BinDesign>& des,
         const std::vector<unsigned>& dset, const std::vector<unsigned>& bin_ids,
         const Eigen::ArrayXd& support);
  
  /// @brief Empty constructor required for technical reasons
  Binner() =default;
  
  /// @brief Apply the binner matrix to a data vector. Returns sum of data in each bin
  Eigen::ArrayXd bin(const Eigen::ArrayXd& data) const;
  
  /// @brief Apply the binner matrix transpose to a bins vector
  Eigen::ArrayXd project_back(const Eigen::ArrayXd& binned) const;
  
  /// @brief Get total number of bins in this binner (even those not used by a data point)
  unsigned get_nbins() const {
    return (data_from_bin_.size()>0) ? (1+data_from_bin_.rbegin()->first) : 0; //bins go from 0 to N
  }
  
  /// @brief Get number of data points in this binner
  unsigned get_ndata() const {
    return (bin_from_data_.size()>0) ? (1+bin_from_data_.rbegin()->first) : 0; 
  }
  
  /// @brief Provide bin id for a bin
  /*! @param binno is the row of the binner matrix, starting at 0. It is specific 
   *                 to a bin id and a dataset.
   *  @return The bin id. It refers to the position within the bin design, which is
   *          independent of the number of datasets.
   */
  unsigned get_bin_id(unsigned binno) const { return bin_ids_[binno]; }
  
  /// @brief Provide dataset id for a bin
  /*! @param binno is the row of the binner matrix, starting at 0. It is specific 
   *                 to a bin id and a dataset.
   *  @return Bins are constructed by binning a dataset. It is this id which is returned.
   */
  unsigned get_bin_dataset(unsigned binno) const { return dset_[binno]; }
  
  /// @brief return sorted list of binner row numbers, each row being non-zero
  std::vector<unsigned> get_used_bins() const;
  
  /// @brief return sorted list of dataset IDs used in this binner
  std::vector<unsigned> get_dataset_ids() const;
  
  /// @brief Print information on binner content, useful for debugging
  void print(int nmax=-1) const;
  
  /// @brief Return the Binner matrix in sparse matrix form
  types::SpMat toSparse() const;
  
  /// @brief extend this binner matrix by adding the second along the diagonal
  void extend_diagonal(const Binner&);
  
  METHLASSO_GET_SET_CONSTREF_DESCR(std::shared_ptr<BinDesign>, bin_design, design used for this binner)
  METHLASSO_GET_SET_OWN_DESCR(std::vector<unsigned>, dset, the dataset ID for each bin)
  METHLASSO_GET_SET_OWN_DESCR(std::vector<unsigned>, bin_ids, the bin IDs)
  typedef std::map<unsigned,unsigned> bin_from_data_t;
  METHLASSO_GET_SET_OWN_DESCR(bin_from_data_t, bin_from_data, a map that associates each data point to a bin)
  typedef std::map<unsigned,std::vector<unsigned> > data_from_bin_t;
  METHLASSO_GET_SET_OWN_DESCR(data_from_bin_t, data_from_bin, a map that associates each bin to its data points)
  
  METHLASSO_GET_SET_OWN_DESCR(Eigen::ArrayXd, support, support is a vector of size Nbins giving the value of the independent variable along these bins)
};
}
}

#endif
