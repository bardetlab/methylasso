#include <algorithm>
#include <Rcpp.h>

#include "Binner.hpp"

namespace MethyLasso {
namespace binned {

Binner::Binner(const types::SpMat& binner, const std::shared_ptr<BinDesign>& des,
               const std::vector<unsigned>& dataset,
               const std::vector<unsigned>& info,
               const Eigen::ArrayXd& support) : bin_design_(des), dset_(dataset), bin_ids_(info), support_(support) {
  std::map<unsigned,unsigned> fwd;
  std::map<unsigned,std::vector<unsigned> > rev;
  unsigned data_per_bin = binner.cols()/double(binner.rows())+1;
  for (int k=0; k < binner.outerSize(); ++k) {
    for (types::SpMat::InnerIterator it(binner,k); it; ++it) {
      if (fwd.emplace(it.col(),it.row()).second == false)
        throw std::invalid_argument("Data point shared across multiple bins!");
      auto& bin = rev[it.row()]; // auto drops reference! Need to put it back...
      if (bin.size()==0) bin.reserve(data_per_bin);
      bin.push_back(it.col());
    }
  }
  std::swap(bin_from_data_, fwd);
  std::swap(data_from_bin_, rev);
}

Eigen::ArrayXd Binner::bin(const Eigen::ArrayXd& data) const {
  unsigned ndata = get_ndata();
  if (ndata != data.size()) throw std::invalid_argument("Invalid size for data vector");
  unsigned nbins = get_nbins();
  Eigen::ArrayXd ret(Eigen::ArrayXd::Zero(nbins));
  for (unsigned i=0; i<ndata; ++i) ret(bin_from_data_.at(i)) += data(i);
  return ret;
}

Eigen::ArrayXd Binner::project_back(const Eigen::ArrayXd& binned) const {
  unsigned nbins = get_nbins();
  if (nbins != binned.size()) throw std::invalid_argument("Invalid size for binned vector");
  unsigned ndata = get_ndata();
  Eigen::ArrayXd ret(Eigen::ArrayXd::Zero(ndata));
  for (auto mit=data_from_bin_.cbegin(); mit != data_from_bin_.cend(); ++mit) {
    auto bin = mit->second;
    for (auto it=bin.cbegin(); it!=bin.cend(); ++it) ret(*it) = binned(mit->first);
  }
  return ret;
}

std::vector<unsigned> Binner::get_used_bins() const {
  std::vector<unsigned> ret;
  ret.reserve(get_nbins());
  for (auto cit = data_from_bin_.cbegin(); cit != data_from_bin_.cend(); ++cit) ret.push_back(cit->first);
  std::sort(ret.begin(),ret.end());
  return ret;
}

std::vector<unsigned> Binner::get_dataset_ids() const {
  std::vector<unsigned> dset(dset_);
  std::sort(dset.begin(), dset.end());
  auto last = std::unique(dset.begin(), dset.end());
  dset.erase(last, dset.end());
  return dset;
}

void Binner::print(int nmax) const {
  unsigned ndata=get_ndata();
  Rcpp::Rcout << "This binner is " << get_nbins() << " x " << ndata << "\n";
  Rcpp::Rcout << "bin_from_data\n";
  if (nmax>0 && ndata > nmax ) ndata = nmax;
  for (unsigned i=0; i < ndata; ++i) Rcpp::Rcout << " data " << i << " in bin " << bin_from_data_.at(i) << "\n";
  Rcpp::Rcout << "data_from_bin\n";
  unsigned reached_data=0;
  for (auto mit=data_from_bin_.cbegin(); mit != data_from_bin_.cend(); ++mit) {
    Rcpp::Rcout << " bin " << mit->first << " with dset " << dset_[mit->first]
                << " and ID " << bin_ids_[mit->first] << " contains data points ";
    auto bin = mit->second;
    for (auto it=bin.cbegin(); it!=bin.cend(); ++it) Rcpp::Rcout << (*it) << " ";
    Rcpp::Rcout << "\n";
    reached_data += bin.size();
    if (nmax>0 && reached_data > nmax) break;
  }
}

types::SpMat Binner::toSparse() const {
  std::vector<Eigen::Triplet<double> > tripletList;
  for (auto cit = bin_from_data_.cbegin(); cit != bin_from_data_.cend(); ++cit)
    tripletList.push_back(Eigen::Triplet<double>(cit->second,cit->first,1.));
  types::SpMat ret(get_nbins(),get_ndata());
  ret.setFromTriplets(tripletList.begin(), tripletList.end());
  return ret;
}


void Binner::extend_diagonal(const Binner& Rhs) {
  if (bin_design_ != Rhs.bin_design_) throw std::invalid_argument("Cannot build block diagonal binner matrix if designs are different");
  //get numbers
  unsigned nbins_left = get_nbins(), nbins_right = Rhs.get_nbins();
  unsigned ndata_left = get_ndata(), ndata_right = Rhs.get_ndata();
  //add Rhs data, shifting indices by the current number of data
  for (auto it = Rhs.bin_from_data_.cbegin(); it != Rhs.bin_from_data_.cend(); ++it)
    bin_from_data_.emplace(ndata_left+it->first, nbins_left+it->second);
  for (auto it = Rhs.data_from_bin_.cbegin(); it != Rhs.data_from_bin_.cend(); ++it) {
    std::vector<unsigned> bin = it->second;
    std::transform(bin.begin(), bin.end(), bin.begin(), [&](unsigned x) -> unsigned { return x+ndata_left; } ); // in place
    data_from_bin_[nbins_left + it->first] = bin;
  }
  //update bin_ids and dset (shifting dset id)
  auto tmp(Rhs.dset_);
  unsigned max_dset(*std::max_element(dset_.begin(),dset_.end()));
  std::transform(tmp.begin(), tmp.end(), tmp.begin(), [&] (unsigned x) -> unsigned { return x+max_dset; } ); //assumes x>0
  dset_.insert(dset_.end(),tmp.begin(),tmp.end());
  bin_ids_.insert(bin_ids_.end(),Rhs.bin_ids_.begin(),Rhs.bin_ids_.end());
  //update support
  Eigen::ArrayXd support(support_.size()+Rhs.support_.size());
  support << support_, Rhs.support_;
  std::swap(support_,support);
}

}
}

