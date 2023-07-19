
template<typename Leg, typename Method>
void Estimator<Leg,Method>::update_binned_residuals(const predictor::IRLSResiduals& z) {
  auto new_binned_residuals = binned::bin_irls_residuals(z,get_binner());
  set_binned_residuals(new_binned_residuals);

  //print debug info
  if (estimator::EstimatorTraits<Leg,Method>::debug_summarizer) {
    Rcpp::Rcout << " estimator " << (this) << " update_binned_residuals\nsizes: binned_residuals.get_zhat()= " << new_binned_residuals.get_zhat().size()
                << " binned_residuals.get_weight()= " << new_binned_residuals.get_weight().size()
                << "\n";
    
    Rcpp::Rcout << "binned_residuals.get_zhat() binned_residuals.get_weight()\n";
    Rcpp::Rcout << (Eigen::MatrixXd(new_binned_residuals.get_zhat().rows(),2)
                << new_binned_residuals.get_zhat(), new_binned_residuals.get_weight()
                ).finished().topRows(std::min<int>(10,new_binned_residuals.get_zhat().rows())) << "\n";
  }
}
  
template<typename Leg, typename Method>
void Estimator<Leg,Method>::update_params() {
  set_old_params(get_params());
  params::Pseudodata pseudo(params::design_pseudodata(get_binned_residuals(), get_design(), get_params()));
  set_pseudodata(pseudo);
  auto hyper = get_hyper();
  auto par = get_fitter().fit_params(pseudo, hyper);
  if (estimator::EstimatorTraits<Leg,Method>::center_params) par = params::center_params(par, pseudo);
  set_params(par);
  
  //print debug info
  if (estimator::EstimatorTraits<Leg,Method>::debug_fitter) {
    Rcpp::Rcout << " estimator " << (this) << " update_params\n sizes: pseudo.get_betahat()= " << pseudo.get_betahat().size()
                << " pseudo.get_pseudoweight()= " << pseudo.get_pseudoweight().size()
                << "\n";
    
    Rcpp::Rcout << "betahat pseudoweight beta\n";
    Rcpp::Rcout << (Eigen::MatrixXd(pseudo.get_betahat().rows(),3)
                << pseudo.get_betahat(), pseudo.get_pseudoweight(), par.get_beta()
                   ).finished().topRows(std::min<int>(10,pseudo.get_betahat().rows())) << "\n";
  }
}

template<typename Leg, typename Method>
void Estimator<Leg,Method>::damp_params(double new_param_weight) {
    auto par = get_fitter().damp_params(get_params(), get_old_params(), new_param_weight);
    set_params(par);
    
    if (estimator::EstimatorTraits<Leg,Method>::debug_fitter) {
      Rcpp::Rcout << "   damp_params(): damping with weight = " << new_param_weight << "\n";
  }
}

template<typename Leg, typename Method>
void Estimator<Leg,Method>::update_hyperparams() {
  if (get_update_hyper()) {
    set_old_hyper(get_hyper());
    auto hyper = get_fitter().fit_hyperparams(get_params());
    set_hyper(hyper);
  }
  
  //print debug info
  if (estimator::EstimatorTraits<Leg,Method>::debug_fitter) {
    Rcpp::Rcout << "   update_hyperparams(): new hyperparameters = " << get_hyper().get_all().transpose() << "\n";
  }
}

template<typename Leg, typename Method>
void Estimator<Leg,Method>::damp_hyperparams(double new_param_weight) {
    auto hyper = get_fitter().damp_hyperparams(get_hyper(), get_old_hyper(), new_param_weight);
    set_hyper(hyper);
    if (estimator::EstimatorTraits<Leg,Method>::debug_fitter) {
      Rcpp::Rcout << "   damp_hyperparams(): damping with weight = " << new_param_weight << "\n";
    }
}

  
