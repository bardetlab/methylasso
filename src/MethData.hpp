#ifndef METH_DATA_HPP
#define METH_DATA_HPP

#include <Rcpp.h>

#include <Eigen/Dense>
#include <vector>

#include "Observations.hpp"

namespace MethyLasso {

class MethData : public data::Observations {
public:
  MethData(const Eigen::ArrayXi& dset, const Eigen::ArrayXi& condition, const Eigen::ArrayXi& replicate,
           const Eigen::ArrayXi& position, const Eigen::ArrayXi& Nmethyl, const Eigen::ArrayXi& Ntot) :
  Observations(dset, condition, replicate, position.cast<double>(), Nmethyl.cast<double>()/Ntot.cast<double>(), Ntot.cast<double>()) {
    if ( (Nmethyl<0).any() || (Ntot<0).any() ) throw std::invalid_argument("Read counts should be >=0 !");
    if ( (Nmethyl>Ntot).any() ) throw std::invalid_argument("Number of methylated reads cannot be larger than number of total reads!");
  }
};





}

#endif

