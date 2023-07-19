#include <cmath>
#include "methylation_helpers.hpp"

namespace MethyLasso {

Rcpp::DataFrame segment_methylation_helper(const Rcpp::DataFrame& df, double tol_val) {
  const Eigen::ArrayXi idx = df["estimate_id"];
  const Eigen::ArrayXi start = df["start"];
  const Eigen::ArrayXd beta = df["beta"];
  const Eigen::ArrayXd betahat = df["betahat"];
  const Eigen::ArrayXd pseudoweight = df["pseudoweight"];
  //
  unsigned curr_idx = idx(0), curr_segment = 0; // current index and segment
  double curr_beta = beta(0); //current value of beta
  double curr_residual = pseudoweight(0) * (betahat(0)-beta(0)); // current residual
  double curr_pseudoweight = pseudoweight(0); //current pseudoweight
  double curr_var = curr_residual * (betahat(0)-beta(0)); //current variance
  std::vector<double> segment_beta, segment_betahat, segment_pseudoweight, segment_var;
  std::vector<unsigned> start_min, start_max, segment_no, segment_idx;
  start_min.push_back(start(0));
  for (unsigned i=1; i<idx.size(); ++i) {
    if (idx(i) == curr_idx) { //we are within the same dataset
      if (std::abs(beta(i)-curr_beta)>tol_val) { //new segment => store old one
        start_min.push_back(start(i));
        start_max.push_back(start(i)-1);
        segment_no.push_back(++curr_segment);
        segment_idx.push_back(curr_idx);
        segment_beta.push_back(curr_beta);
        segment_pseudoweight.push_back(curr_pseudoweight);
        segment_betahat.push_back(curr_residual/curr_pseudoweight + curr_beta);
        segment_var.push_back(curr_var/curr_pseudoweight-curr_residual*curr_residual/(curr_pseudoweight*curr_pseudoweight));
        //
        curr_beta = beta(i);
        curr_residual = pseudoweight(i) * (betahat(i)-beta(i));
        curr_var = curr_residual * (betahat(i)-beta(i));
        curr_pseudoweight = pseudoweight(i);
      } else { //accumulate residuals and pseudoweight
        curr_residual += pseudoweight(i) * (betahat(i)-beta(i));
        curr_var += pseudoweight(i) * (betahat(i)-beta(i)) * (betahat(i)-beta(i));
        curr_pseudoweight += pseudoweight(i);
      }
    } else { //dataset change
      start_min.push_back(start(i));
      start_max.push_back(start(i-1));
      segment_no.push_back(++curr_segment);
      segment_idx.push_back(curr_idx);
      segment_beta.push_back(curr_beta);
      segment_pseudoweight.push_back(curr_pseudoweight);
      segment_betahat.push_back(curr_residual/curr_pseudoweight + curr_beta);
      segment_var.push_back(curr_var/curr_pseudoweight-curr_residual*curr_residual/(curr_pseudoweight*curr_pseudoweight));
      //
      curr_idx = idx(i);
      curr_segment = 0;
      curr_beta = beta(i);
      curr_residual = pseudoweight(i) * (betahat(i)-beta(i));
      curr_var = curr_residual * (betahat(i)-beta(i));
      curr_pseudoweight = pseudoweight(i);
    }
  }
  start_max.push_back(start(idx.size()-1));
  segment_no.push_back(++curr_segment);
  segment_idx.push_back(curr_idx);
  segment_beta.push_back(curr_beta);
  segment_pseudoweight.push_back(curr_pseudoweight);
  segment_betahat.push_back(curr_residual/curr_pseudoweight + curr_beta);
  segment_var.push_back(curr_var/curr_pseudoweight-curr_residual*curr_residual/(curr_pseudoweight*curr_pseudoweight));
  for (unsigned i=0; i<segment_var.size(); ++i) segment_var[i] = std::sqrt(segment_var[i]);
  return Rcpp::DataFrame::create(_["estimate_id"]=segment_idx, _["segment_id"]=segment_no,
                                 _["start"]=start_min, _["end"]=start_max, _["beta"]=segment_beta,
                                 _["betahat"]=segment_betahat, _["pseudoweight"]=segment_pseudoweight,
                                 _["stdev"]=segment_var);
}

}


