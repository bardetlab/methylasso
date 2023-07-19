#ifndef BIN_DESIGN_HPP
#define BIN_DESIGN_HPP

#include <Rcpp.h>

/*! \file BinDesign.hpp
 * @brief implementation of BinDesign and friend functions
 */

#include "typedefs.hpp"
#include "macros.hpp"

namespace MethyLasso {
namespace binned {

/// @brief BinDesign contains the boundaries for a given binning
struct BinDesign {
  BinDesign(unsigned Nbins_, const std::map<double,unsigned>& map_, const std::vector<double>& starts_,
            const std::vector<double>& ends_) : Nbins(Nbins_), map(map_), starts(starts_), ends(ends_) {}
  /// number of bins, numbered from 0 to Nbins
  unsigned Nbins;
  /// std::map of bin starts to bin numbers
  std::map<double,unsigned> map;
  /// vectors of bin starts and ends
  std::vector<double> starts, ends;
};

}
}

#endif
