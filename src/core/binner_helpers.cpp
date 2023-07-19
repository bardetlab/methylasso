#include <numeric>
#include <stdexcept>
#include <Rcpp.h>

#include "binner_helpers.hpp"

namespace MethyLasso {
namespace binned {


std::shared_ptr<BinDesign> make_even_bins(double lower, double upper, unsigned Nbins) {
  double size = upper - lower;
  std::map<double,unsigned> ret;
  std::vector<double> starts, ends;
  for (unsigned i=0; i<=Nbins; ++i) { //size of map is Nbins+1
    double start = lower + i*size/double(Nbins);
    ret[start]=i;
    if (i<Nbins) starts.push_back(start);
    if (i>0) ends.push_back(start);
  }
  return std::make_shared<BinDesign>(Nbins,ret,starts,ends);
}

std::shared_ptr<BinDesign> make_unique_bins(const Eigen::ArrayXi& pos) {
  //get unique values
  std::vector<unsigned> unique_pos(pos.data(), pos.data() + pos.rows());
  std::sort(unique_pos.begin(), unique_pos.end());
  auto unique_it = std::unique(unique_pos.begin(), unique_pos.end());
  unique_pos.resize( std::distance(unique_pos.begin(), unique_it) );
  //build BinDesign
  const unsigned Nbins = unique_pos.size();
  std::map<double,unsigned> ret;
  std::vector<double> starts, ends;
  for (unsigned i=0; i < Nbins; ++i) {
    double dpos = unique_pos[i];
    ret[dpos]=i;
    if (i<Nbins) starts.push_back(dpos);
    if (i>0) ends.push_back(dpos);
  }
  ret[unique_pos.back()+1] = Nbins; //map is size N+1
  ends.push_back(unique_pos.back()+1);
  return std::make_shared<BinDesign>(Nbins,ret,starts,ends);
}



std::shared_ptr<binned::Binner> make_binner_using_bin_design(const Eigen::ArrayXd& data, const std::shared_ptr<BinDesign>& design, bool drop) {
  //first construct the full triplet list
  unsigned Ndata = data.rows();
  std::vector<Eigen::Triplet<double> > tripletList;
  tripletList.reserve(Ndata);
  for (unsigned i=0; i<Ndata; ++i) {
    auto it = design->map.upper_bound(data(i));
    unsigned j;
    if (it == design->map.end()) {
      j = design->map.rbegin()->second -1; //place in last bin
    } else {
      j = it->second - 1; //upper_bound returns pointer to next bin
    }
    tripletList.push_back(Eigen::Triplet<double>(j,i,1.));
  }
  //prepare bin information
  std::vector<unsigned> bin_ids(design->Nbins);
  std::iota(std::begin(bin_ids),std::end(bin_ids),0); //set bin_ids to 0, 1, ..., Nbins-1
  //remove unused bins if required
  unsigned final_num_bins = design->Nbins;
  if (drop) {
    //get which rows are empty
    std::vector<bool> has_value(design->Nbins,false);
    for (auto tr : tripletList) has_value[tr.row()] = true; 
    //adapt bin info
    std::vector<unsigned> new_bin_ids;
    for (unsigned old_idx=0; old_idx<design->Nbins; old_idx++) if(has_value[old_idx]) new_bin_ids.push_back(old_idx);
    swap(new_bin_ids,bin_ids);
    //create map from old to new indices
    std::map<unsigned,unsigned> row_map;
    unsigned new_idx=0;
    for (unsigned old_idx=0; old_idx<design->Nbins; old_idx++) if(has_value[old_idx]) row_map[old_idx]=new_idx++;
    //make new triplet list, dropping empty rows
    std::vector<Eigen::Triplet<double> > newTripletList;
    newTripletList.reserve(Ndata);
    for (auto tr : tripletList) newTripletList.push_back(Eigen::Triplet<double>(row_map[tr.row()],tr.col(),tr.value()));
    swap(tripletList,newTripletList);
    final_num_bins = new_idx;
  }
  //compute support from bin starts
  Eigen::Map<const Eigen::ArrayXd> support(design->starts.data(), design->starts.size());
  //Rcpp::Rcout << "Binner support (sz= " << support.size() << " ) : " << support.head(5).transpose() << " ... " << support.tail(5).transpose() << "\n";
  //now construct binner
  std::vector<unsigned> dset(final_num_bins,1); //set dataset id to 1 (cannot be 0 otherwise breaks block_diagonal)
  types::SpMat mat(final_num_bins,Ndata);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  mat.makeCompressed();
  return std::make_shared<Binner>(mat, design, dset, bin_ids, support);
}

std::vector<Eigen::ArrayXd> split_by_dataset(const Eigen::ArrayXi& dset, const Eigen::ArrayXd& pos) {
  //make list of dataset starts and sizes
  const unsigned npos = pos.size();
  std::vector<unsigned> dset_starts(1,0),dset_sizes;
  for (unsigned i=0; i<npos; ++i) {
    if (dset(i) != dset_starts.size()) {
      dset_sizes.push_back(i-dset_starts.back());
      dset_starts.push_back(i);
    }
  }
  dset_sizes.push_back(npos-dset_starts.back());
  //split data vector
  std::vector<Eigen::ArrayXd> ret;
  for (unsigned i=0; i<dset_starts.size(); ++i) ret.push_back(pos.segment(dset_starts[i],dset_sizes[i]));
  return ret;
}

}
}

