#include <Eigen/Sparse>
#include <stdexcept>
#include <Rcpp.h>

#include "design_helpers.hpp"
#include "typedefs.hpp"
#include "Binner.hpp"
#include "Design.hpp"
#include "DataDesign.hpp"

namespace MethyLasso {
namespace params {

Design make_identity_design(const binned::Binner& binner) {
  unsigned nbins = binner.get_nbins();
  auto X = types::SpMat(nbins,nbins);
  X.setIdentity();
  auto starts = binner.get_bin_design()->starts; //must store intermediate
  Eigen::Map<const Eigen::ArrayXd> support(starts.data(), starts.size());
  return Design(X,support);
}

Design make_group_design(const binned::Binner& binner, const data::DataDesign& grouped_design) {
  //get design parameters
  auto bin_nos = binner.get_used_bins(); //skips empty rows
  unsigned nbins = binner.get_nbins(); //includes empty rows (if modelled in the binner)
  unsigned ntotal = binner.get_bin_design()->Nbins; //number of parameters per dataset
  Eigen::SparseVector<double> coefs(grouped_design.get_X()); // vector of coefficients
  std::vector<bool> is_nonzero(coefs.size(),false); //whether coef is zero or not
  for (Eigen::SparseVector<double>::InnerIterator it(coefs); it; ++it) is_nonzero[it.index()]=true;
  //then build triplets
  std::vector<Eigen::Triplet<double> > tripletList;
  for (auto cit = bin_nos.cbegin(); cit != bin_nos.cend(); ++cit) {
    //get bin info
    unsigned bin_id = binner.get_bin_id(*cit);
    unsigned bin_dset = binner.get_bin_dataset(*cit);
    //write coef for that dataset
    if (is_nonzero[bin_dset-1]) tripletList.push_back(Eigen::Triplet<double>(*cit,bin_id,coefs.coeffRef(bin_dset-1)));
  }
  //and fill the matrix
  Eigen::SparseMatrix<double> ret(nbins,ntotal);
  ret.setFromTriplets(tripletList.begin(), tripletList.end());
  //finally, compute support from binner design
  auto starts = binner.get_bin_design()->starts; //must store intermediate
  Eigen::Map<const Eigen::ArrayXd> support(starts.data(), starts.size());
  //Rcpp::Rcout << "Design support (sz= " << support.size() << " ) : " << support.head(5).transpose() << " ... " << support.tail(5).transpose() << "\n";
  return Design(ret,support);
}


}
}
