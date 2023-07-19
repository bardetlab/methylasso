#ifndef OBSERVATIONS_HPP
#define OBSERVATIONS_HPP

#include <Rcpp.h>

#include <Eigen/Dense>
#include <vector>

#include "macros.hpp"

namespace MethyLasso {
namespace data {

/// @brief Observations contains the data which is to be fitted
/*!
 *  To implement your own data container, simply derive from this class. Derived classes
 *  should not contain any data members, since they will be discarded (class is internally cast to Observations)
 */
class Observations {
public:
  Observations(const Eigen::ArrayXi& dataset, const Eigen::ArrayXi& condition, const Eigen::ArrayXi& replicate,
               const Eigen::ArrayXd& position, const Eigen::ArrayXd& observed, const Eigen::ArrayXd& weight) :
  dataset_(dataset), condition_(condition), replicate_(replicate), position_(position), observed_(observed), Nobs_(weight), Ndata_(position_.size()) {
      check_observations();
  }
  
private:

  void check_observations() const {
      //number of datapoints >= 2
      if (Ndata_<=1) throw std::invalid_argument("Need at least 2 data points!");
      //all arrays of same size
      if (dataset_.size()!=Ndata_ || condition_.size()!=Ndata_ || replicate_.size()!=Ndata_ || position_.size()!=Ndata_ || observed_.size()!=Ndata_ || Nobs_.size()!=Ndata_)
          throw std::invalid_argument("All input arrays must be of same size");
      double last_pos = position_(0);
      unsigned last_dset = dataset_(0);
      unsigned last_condition = condition_(0);
      unsigned last_replicate = replicate_(0);
      //dataset starts at 1
      if (last_dset != 1) throw std::invalid_argument("First dataset ID must be 1");
      if (last_condition != 1) throw std::invalid_argument("First condition ID must be 1");
      if (last_replicate != 1) throw std::invalid_argument("First replicate ID must be 1");
      //data must be sorted by dataset ID and position
      for (unsigned i=1; i<Ndata_; ++i) {
          if (last_dset == dataset_(i)) {
              if (last_pos >= position_(i)) throw std::invalid_argument("data is not sorted by dataset and position");
          } else {
              if (dataset_(i) != ++last_dset) throw std::invalid_argument("dataset IDs must start at 1 and increase without gaps");
              if (condition_(i) == last_condition) {
                if (replicate_(i) != ++last_replicate) throw std::invalid_argument("replicate ID must increase without gaps within the same condition");
              } else {
                if (condition_(i) != ++last_condition) throw std::invalid_argument("condition ID must start at 1 and increase without gaps");
                last_replicate = replicate_(i);
                if (last_replicate != 1) throw std::invalid_argument("replicate ID must start at 1 within the same condition");
              }
          }
          last_pos = position_(i);
      }
  }
  
  //these variables must be implemented for Posterior to work
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXi, dataset, dataset index)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXi, condition, condition index)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXi, replicate, replicate index)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXd, position, independent variable)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXd, observed, dependent variable)
  METHLASSO_GET_CONST_OWN_DESCR(Eigen::ArrayXd, Nobs, weight of each dependent variable or number of draws)
  METHLASSO_GET_CONST_SIMPLE_DESCR(unsigned, Ndata, number of datapoints)
  
};


}
}

#endif

