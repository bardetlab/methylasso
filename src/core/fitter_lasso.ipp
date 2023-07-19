template<typename Library>
params::Params<Lasso> Fitter<Lasso,Library>::fit_params(const params::Pseudodata& pseudo, const params::Hyperparams<Lasso>& hyper) const {
  //we do not optimize lambda2 here
  if (hyper.get_lambda2()<=0) throw std::invalid_argument("lambda2 should be >0");
  
  //update beta
  Eigen::ArrayXd betahat = pseudo.get_betahat(); //essential, because the std::vector uses that memory location
  Eigen::ArrayXd pseudoweight = pseudo.get_pseudoweight();
  std::vector<double> y(betahat.data(), betahat.data() + betahat.size());  
  std::vector<double> wt(pseudoweight.data(), pseudoweight.data() + pseudoweight.size());  
  flo_.optimize(y,wt,hyper.get_lambda2());
  std::vector<double> beta_std(flo_.get(0, get_lambda1()));
  Eigen::Map<Eigen::ArrayXd> beta(beta_std.data(), beta_std.size());
  return params::Params<Lasso>{beta};
}

template<typename Library>
params::Hyperparams<Lasso> Fitter<Lasso,Library>::fit_hyperparams(const params::Params<Lasso>& par) const {
  Eigen::ArrayXd beta(par.get_beta());
  std::vector<double> beta_std(beta.data(), beta.data() + beta.size());  
  double lambda2 = flo_.optimize_lambda(beta_std);
  return params::Hyperparams<Lasso>{lambda2};
}

 
/// @brief return the sum of minus log prior and hyperprior terms
template<typename Library>
double Fitter<Lasso,Library>::get_minus_log_prior(const params::Params<Lasso>& par,
                                              const params::Hyperparams<Lasso>& hyper) const {
  Eigen::ArrayXd beta(par.get_beta());
  std::vector<double> beta_std(beta.data(), beta.data() + beta.size());  
  double lambda2 = hyper.get_lambda2();
  return flo_.get_minus_log_prior(beta_std,lambda2);
}


