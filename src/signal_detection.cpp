#include "signal_detection.hpp"
#include "NormalDistribution.hpp"
#include "BetaBinomialRateDistribution.hpp"
#include "link_functions.hpp"
#include "binned_genome.hpp"
#include "unique_positions.hpp"
#include "methylation_detection.hpp"

namespace MethyLasso {

Rcpp::List signal_detection_fast_R(const Rcpp::DataFrame& data, double lambda2,
                                   double tol_val, unsigned max_iter) {
  if (! (data.containsElementNamed("dataset") && data.containsElementNamed("condition") && data.containsElementNamed("replicate")
           && data.containsElementNamed("pos") && data.containsElementNamed("Nmeth") && data.containsElementNamed("coverage")) )
    Rcpp::stop("data should contain columns dataset, condition, replicate, pos, Nmeth and coverage!");
  //
  // normal distribution with one basis function per CpG (distance-independent)
  //
  MethData obs(data["dataset"],data["condition"],data["replicate"],data["pos"],data["Nmeth"],data["coverage"]);
  //
  const double max_lambda2 = lambda2 + 1;
  const double lambda1=0;
  const unsigned fake_nouter = 1;
  UniquePositionsConfig mconf(lambda2, max_lambda2, lambda1, max_iter, tol_val, fake_nouter);
  //
  const double fake_hyper=1;
  const bool optimize_hyper=false;
  const bool optimize_lambda2=false;
  //
  const unsigned num_dsets = obs.get_dataset().maxCoeff();
  const unsigned num_conditions = obs.get_condition().maxCoeff();
  //Rcpp::Rcout << "Fast signal detection with " << num_dsets << " datasets and " << num_conditions << " conditions"
              //<< " using lambda2 = " << lambda2 << "\nDesign matrix:\n";
  auto design = build_signal_design(obs);
  //Rcpp::Rcout << design.get_X() << "\n";
  //
  auto ret = methylation_detection<distribution::NormalDistribution,
                                   link::Identity, UniquePositions,
                                   AdaptiveDampingIRLSPolicy<1,5>, BasicPriorPolicy, BasicHyperPolicy
    >(obs, design, mconf, optimize_lambda2, fake_hyper, optimize_hyper);
  //
  //Rcpp::Rcout << "The optimization has " << (Rcpp::as<bool>(ret["converged"])?"":"NOT ") << "converged\n";
  //Rcpp::Rcout << "Done.\n";
  return Rcpp::List::create(_["mu"]=ret["mu"], _["lambda2"]=ret["lambda2"], _["converged"]=ret["converged"],
                            _["offset"]=ret["offset"], _["estimates"]=ret["estimates"],
                            _["detection.type"]="signal", _["model"]="normal");
}

Rcpp::List signal_detection_R(const Rcpp::DataFrame& data, bool optimize_lambda2, double lambda2, double max_lambda2,
                              bool optimize_hyper, double hyper, unsigned nouter, double bf_per_kb, double tol_val,
                              unsigned max_iter) {
  if (! (data.containsElementNamed("dataset") && data.containsElementNamed("condition") && data.containsElementNamed("replicate")
           && data.containsElementNamed("pos") && data.containsElementNamed("Nmeth") && data.containsElementNamed("coverage")) )
    Rcpp::stop("data should contain columns dataset, condition, replicate, pos, Nmeth and coverage!");
  // beta binomial with CpG-distance-dependent modelling
  MethData obs(data["dataset"],data["condition"],data["replicate"],data["pos"],data["Nmeth"],data["coverage"]);
  //
  double lambda1=0;
  BinnedGenomeConfig mconf(lambda2, max_lambda2, lambda1, bf_per_kb, max_iter, tol_val, nouter);
  //
  const unsigned num_dsets = obs.get_dataset().maxCoeff();
  const unsigned num_conditions = obs.get_condition().maxCoeff();
  //Rcpp::Rcout << "Signal detection with " << num_dsets << " datasets and " << num_conditions << " conditions"
  //            << " using lambda2 = " << lambda2 << "\nDesign matrix:\n";
  auto design = build_signal_design(obs);
  //Rcpp::Rcout << design.get_X() << "\n";
  //
  auto ret = methylation_detection<distribution::BetaBinomialRateDistribution,
                                   link::Logit, BinnedGenome,
                                   AdaptiveDampingIRLSPolicy<1,5>,
                                   AdaptiveDampingPriorPolicy<1,5>,
                                   AdaptiveDampingHyperPolicy<1,5>
                         >(obs, design, mconf, optimize_lambda2, hyper, optimize_hyper);
  //
  //Rcpp::Rcout << "The optimization has " << (Rcpp::as<bool>(ret["converged"])?"":"NOT ") << "converged\n";
  //Rcpp::Rcout << "Done.\n";
  return Rcpp::List::create(_["mu"]=ret["mu"], _["lambda2"]=ret["lambda2"], _["hyper"]=ret["hyper"],
                            _["converged"]=ret["converged"], _["offset"]=ret["offset"], _["estimates"]=ret["estimates"],
                            _["detection.type"]="signal", _["model"]="beta-binomial");
  
}

}


