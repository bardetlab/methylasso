template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class... Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::initial_guess(const data::Observations& data, Est&&... est) {
  predictor::MeanVector mu0 = distribution_.get_mean_initial_guess(data);
  mu0 = link::restrict_to_support<Link>(mu0);
  predictor::Moments moments(predictor::get_moments_from_distribution(mu0,data,distribution_));
  /*Rcpp::Rcout << "obs: " << data.get_obs().transpose().segment(213,3) << "\n";
  Rcpp::Rcout << "Nobs: " << data.get_Nobs().transpose().segment(213,3) << "\n";
  Rcpp::Rcout << "mu0: " << moments.get_mean().transpose().segment(213,3) << "\n";
  Rcpp::Rcout << "V0: " << moments.get_variance().transpose().segment(213,3) << "\n";*/
  predictor::IRLSResiduals z = link::transform_initial_guess<Link>(moments);
  /*Rcpp::Rcout << "z.residuals: " << z.residuals.transpose().segment(213,3) << "\n";
  Rcpp::Rcout << "z.weights: " << z.weights.transpose().segment(213,3) << "\n";*/
  update_all_binned_residuals(z, std::forward<Est>(est)...);
  set_converged(false);
  //check for user interrupt
  Rcpp::checkUserInterrupt();
}
  
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class... Est>
predictor::MeanVector Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::iterate(const data::Observations& data, Est&&... est) {
  predictor::MeanVector mu,mu_old;
  set_converged(false);
  double old_mlogp = get_minus_log_posterior(data,std::forward<Est>(est)...)(0);
  double mlogp;
  //Rcpp::Rcout << "   -log(posterior) init is " << old_mlogp << "\n";
  for (unsigned i=0; i<get_maxiter(); ++i) {
    //update and adapt params
    update_all_params(std::forward<Est>(est)...);
    adapt_all_params(i, old_mlogp, data, std::forward<Est>(est)...);
    //check for convergence
    mu = get_mu(std::forward<Est>(est)...);
    mlogp = get_minus_log_posterior(data,std::forward<Est>(est)...)(0);
    //Rcpp::Rcout << "   -log(posterior) step " << (i+1) << " is " << mlogp << "\n";
    if (std::isnan(mlogp)) Rcpp::stop("NaN found, aborting...");
    if (i>0) { //otherwise mu_old is undefined
      double diff = predictor::get_max_absolute_deviation(mu,mu_old);
      //Rcpp::Rcout << "   Iteration " << (i+1) << " reached tolerance " << diff << " (criterion : " << get_tol_val() << ")\n";
      if (diff < get_tol_val()) {
        set_converged(true);
        //Rcpp::Rcout << "   IRLS converged after " << (i+1) << " / " << maxiter_ << " steps\n";
        break;
      }
    }
    mu_old = mu;
    old_mlogp = mlogp;
    //update summaries
    predictor::Moments moments(predictor::get_moments_from_distribution(mu,data,distribution_));
    predictor::IRLSResiduals z = link::get_irls_residuals<Link>(data, moments);
    update_all_binned_residuals(z, std::forward<Est>(est)...);
    //check for user interrupt
    Rcpp::checkUserInterrupt();
  }
  if ( get_maxiter()>1 && (!get_converged()) ) 
     Rcpp::Rcout << "Warning: did not converge after " << maxiter_ << " steps\n";
  if ( get_maxiter()==1 ) set_converged(true);
  return mu;
}
  
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class... Est>
Eigen::Array3d Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_minus_log_posterior(const data::Observations& data, Est&&... est) const {
  predictor::MeanVector mu(get_mu(std::forward<Est>(est)...));
  double lik = predictor::get_minus_log_pdf(data,mu.get_mu(),distribution_);
  double pri = get_minus_log_prior(std::forward<Est>(est)...);
  Eigen::Array3d ret(lik+pri,lik,pri);
  return ret;
}


template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class... T>
predictor::MeanVector Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_mu(T&&... t) const {
  auto eta = get_eta(std::forward<T>(t)...);
  //Rcpp::Rcout << "eta(12717)= " << eta.get_eta()(12717) << "\n";
  predictor::MeanVector mu = link::get_mu<Link>(eta);
  // restrict to support to avoid problems caused by boundaries
  mu = link::restrict_to_support<Link>(mu);
  //Rcpp::Rcout << "mu(12717)= " << mu.get_mu()(12717) << "\n";
  return mu;
}
  
//parameter update: call update_params() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> 
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_all_params(Est& est, T&&... t) {
  update_all_params(est);
  update_all_params(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_all_params(Est& est) {
  est.update_params();
}

template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Data, class... Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::adapt_all_params(unsigned step, double old_mlogp, const Data& data, Est&&... est) {
  // damp iteratively until the -log(p) decreases again
  unsigned i=0;
  double mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  //Rcpp::Rcout << "param step " << step << " damp " << -1 << " damp target old_mlogp= " << old_mlogp <<  "\n";
  while(IRLSPolicy::condition_check(step,i++,mlogp,old_mlogp)) {
    //Rcpp::Rcout << "param step " << step << " damp " << i << " damped to mlogp= " << mlogp << "\n";
    damp_all_params(std::forward<Est>(est)...);
    mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  }
  //Rcpp::Rcout << "param step " << step << " damp " << i << " damped to mlogp= " << mlogp << "\n";
}

//parameter update: call damp_params() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> 
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::damp_all_params(Est& est, T&&... t) {
  damp_all_params(est);
  damp_all_params(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::damp_all_params(Est& est) {
  est.damp_params(IRLSPolicy::get_new_param_weight());
}
  
//hyperparameter update: call update_hyperparams() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> 
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_prior_hyperparams(Est& est, T&&... t) {
  update_prior_hyperparams(est);
  update_prior_hyperparams(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_prior_hyperparams(Est& est) {
  est.update_hyperparams();
}

template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Data, class... Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::adapt_prior_hyperparams(unsigned step, double old_mlogp, const Data& data, Est&&... est) {
  // damp iteratively until the -log(p) decreases again
  unsigned i=0;
  double mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  Rcpp::Rcout << "hyper step " << step << " damp " << -1 << " damp target old_mlogp= " << old_mlogp <<  "\n";
  while(PriorPolicy::condition_check(step,i++,mlogp,old_mlogp)) {
    Rcpp::Rcout << "hyper step " << step << " damp " << i << " damped to mlogp= " << mlogp << "\n";
    damp_prior_hyperparams(std::forward<Est>(est)...);
    mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  }
  Rcpp::Rcout << "hyper step " << step << " damp " << i << " damped to mlogp= " << mlogp << "\n";
}

//hyperparameter update: call damp_hyperparams() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> 
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::damp_prior_hyperparams(Est& est, T&&... t) {
  damp_prior_hyperparams(est);
  damp_prior_hyperparams(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::damp_prior_hyperparams(Est& est) {
  double param_weight = PriorPolicy::get_new_param_weight();
  est.damp_hyperparams(param_weight);
}

template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Data, class... Est>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::adapt_distribution_hyperparams(unsigned step, double old_mlogp, const Data& data, Est&&... est) {
  // damp iteratively until the -log(p) decreases again
  unsigned i=0;
  double mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  while(HyperPolicy::condition_check(step,i++,mlogp,old_mlogp)) {
    damp_distribution_hyperparams();
    mlogp = get_minus_log_posterior(data, std::forward<Est>(est)...)(0);
  }
}

//hyperparameter update: call damp_hyperparams() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::damp_distribution_hyperparams() {
  double param_weight = HyperPolicy::get_new_param_weight();
  distribution_.damp_hyperparams(param_weight);
}

//get_prior_hyperparams: return vector of all get_hyperparams() calls
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T>
Eigen::ArrayXd Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_prior_hyperparams(const Est& est, T&&... t) const {
  Eigen::ArrayXd a(get_prior_hyperparams(est));
  Eigen::ArrayXd b(get_prior_hyperparams(std::forward<T>(t)...));
  Eigen::ArrayXd ret(a.rows()+b.rows());
  ret << a,b;
  return ret;
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
Eigen::ArrayXd Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_prior_hyperparams(const Est& est) const {
  return est.get_hyperparams();
}

template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class... Est>
predictor::MeanVector Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::maximize_posterior(unsigned nouter, bool optimize_prior,
                                                    bool optimize_hyper, const data::Observations& data, Est&&... est) {
  predictor::MeanVector mu;
  if (nouter>1) {
    predictor::MeanVector mu_old;
    for (unsigned i=0; i<nouter; ++i) {
        //optimize mu
        Rcpp::Rcout << "  Outer iteration " << (i+1) << "\n";
        mu = iterate(data, std::forward<Est>(est)...);
        //check for convergence
        if (i>0) { //otherwise mu_old is undefined
            double diff = predictor::get_max_absolute_deviation(mu,mu_old);
            Rcpp::Rcout << "   Outer iteration " << (i+1) << " reached tolerance " << diff << " (criterion : " << get_tol_val() << ")\n";
            if (diff < get_tol_val()) {
              set_converged(get_converged() && true); //get_converged() includes IRLS convergence info
              Rcpp::Rcout << "   outer iteration converged after " << (i+1) << " / " << nouter << " steps\n";
              break;
            } else {
                set_converged(false);
            }
        }
        mu_old = mu;
        //update hyperparameters
        //betware the order: prior updates mu, which is used by distribution
        double mlogp = get_minus_log_posterior(data,std::forward<Est>(est)...)(0);
        if (optimize_prior) {
            update_prior_hyperparams(std::forward<Est>(est)...);
            adapt_prior_hyperparams(i, mlogp, data, std::forward<Est>(est)...);
            if (optimize_hyper) mlogp = get_minus_log_posterior(data,std::forward<Est>(est)...)(0);
        }
        if (optimize_hyper) {
          distribution_.update_hyperparams(data.get_observed(), data.get_Nobs(),
                                           get_mu(std::forward<Est>(est)...).get_mu());
          adapt_distribution_hyperparams(i, mlogp, data, std::forward<Est>(est)...);
        }
        //check for weird estimates
        bool prior_ok = get_prior_hyperparams(std::forward<Est>(est)...).isFinite().all();
        bool dist_ok = get_distribution_hyperparams().isFinite().all();
        if (!(prior_ok && dist_ok)) {
            set_converged(false);
            Rcpp::Rcout << "   outer iteration leads to infinite/nan values, aborting...\n";
            break;
        }
    }
  } else {
    mu = iterate(data, std::forward<Est>(est)...);
  }
  return mu;
}

//binned_residuals update: call update_summaries(z) for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T>
void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_all_binned_residuals(const predictor::IRLSResiduals& z, Est& est, T&&... t) {
  update_all_binned_residuals(z, est);
  update_all_binned_residuals(z, std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est> void Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::update_all_binned_residuals(const predictor::IRLSResiduals& z, Est& est) {
  est.update_binned_residuals(z);
}
  
//get_eta: return sum of all get_linear_response() calls
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T>
predictor::LinearResponse Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_eta(const Est& est, T&&... t) const {
  return get_eta(est) + get_eta(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est>
predictor::LinearResponse Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_eta(const Est& est) const {
  auto ret = est.get_linear_response();
  //Rcpp::Rcout << "get_linear_response: " << ret.get_eta().tail(3).transpose() << "\n";
  return ret;
}
  
//parameter update: call get_minus_log_prior() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> double Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_minus_log_prior(const Est& est, T&&... t) const {
  return get_minus_log_prior(est) + get_minus_log_prior(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est> double Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_minus_log_prior(const Est& est) const {
  return est.get_minus_log_prior();
}

//parameter update: call get_irls_target_with_prior() for each estimator
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est, class... T> double Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_irls_target_with_prior(const Est& est, T&&... t) const {
  return get_irls_target_with_prior(est) + get_irls_target_with_prior(std::forward<T>(t)...);
}
template<typename Distribution, typename Link, typename IRLSPolicy, typename PriorPolicy, typename HyperPolicy>
template<class Est> double Posterior<Distribution,Link,IRLSPolicy,PriorPolicy,HyperPolicy>::get_irls_target_with_prior(const Est& est) const {
  return est.get_irls_target_with_prior();
}

