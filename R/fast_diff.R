#' @include dependencies.R
NULL


#' Methylation difference detection using fused lasso, for 2 datasets (each can
#' have replicates)
#' 
#' Fast model: normal distribution with identity link and distance-independent modelling
#' 
#' single chromosome
#' 
difference_detection_fast_singlechr = function(data, ref, lambda2=50,
                                               tol_val=0.01, max_iter=1) {
  ret = MethyLasso:::signal_detection_fast_singlechr(data,lambda2,tol_val,max_iter)
  est = as.data.table(ret$estimates)
  est_ref = est[estimate_id==1]
  est = merge(est[estimate_id!=1],
              est_ref[,.(start, beta_ref = est_ref$beta, betahat_ref = betahat,
                         pseudoweight_ref = pseudoweight)], by = "start", all.x=T)
  est = est[, .(estimate_id, start, beta = beta - beta_ref,
                betahat = (betahat * pseudoweight - betahat_ref * pseudoweight_ref) / (pseudoweight + pseudoweight_ref),
                pseudoweight = pseudoweight + pseudoweight_ref)]
  est[pseudoweight==0, betahat:=0]
  est = rbind(est_ref, est)
  setkey(est, estimate_id, start)
  ret$estimates = as.data.frame(est)
  ret$detection.type = "difference"
  ret$ref.cond = ref
  return(ret)
}

