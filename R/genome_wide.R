#' @include dependencies.R
NULL

#' Helper to combine two return values from signal/difference calls
#' 
combine_ret = function(ret1, ret2) {
  stopifnot(ret1$model == ret2$model && ret1$detection.type == ret2$detection.type)
  if (ret1$detection.type== "difference") stopifnot(ret1$ref.cond == ret2$ref.cond)
  ret =  list(mu=c(ret1$mu,ret2$mu), lambda2=c(ret1$lambda2,ret2$lambda2), hyper=c(ret1$hyper,ret2$hyper),
              converged=c(ret1$converged,ret2$converged), offset=c(ret1$offset,ret2$offset),
              estimates = rbind(ret1$estimates, ret2$estimates), chr = c(ret1$chr, ret2$chr),
              detection.type=ret1$detection.type, model = ret1$model)
  if (ret1$detection.type== "difference") ret$ref.cond = ret1$ref.cond
  return(ret)
}

#' Methylation detection using fused lasso, for 1 dataset or several replicates
#'
#' Uses beta binomial distribution with logit link
#' 
#' @param data dataframe containing columns dataset, condition, replicate, chr, pos, Nmethyl and coverage.
#'  dataset, condition and replicate are integer vectors or factors, starting at 1 and increasing without gaps.
#'  chr is the chromosome, as a character or factor. Each chromosome is treated separately.
#'  pos is the genomic coordinate of each CpG, sorted increasingly. Nmethyl is the number
#'  of methylated reads. Ntot is the number of unmethylated reads.
#' @param optimize_lambda2 whether to optimize fusion penalty or not (default FALSE)
#' @param lambda2 the initial value of lambda2 (default 25)
#' @param max_lambda2 the maximum value of lambda2 (default 1000)
#' @param optimize_hyper whether to optimize it (default FALSE)
#' @param hyper the value of the hyperparameter (default 100)
#' @param nouter the number of iterations to optimize hyper and/or lambda2 (default 20)
#' @param bf_per_kb Number of basis functions per kb (default 50)
#' @param tol_val Target precision for the estimated methylation value (default 0.01)
#' @param max_iter maximum number of IRLS iterations for the beta binomial model
#' @param ncores number of cores to parallelize on (default 1)
#' 
#' @export
signal_detection = function(data, optimize_lambda2=F, lambda2=25, max_lambda2=1000,
                            optimize_hyper=F, hyper=100, nouter=20, bf_per_kb=50, tol_val=0.01,
                            max_iter=30, ncores=1) {
  chrs=data[,unique(chr)]
  registerDoParallel(cores=ncores)
  ret = foreach (ch=chrs,.combine=MethyLasso:::combine_ret, .errorhandling = "remove") %dopar% {
    mdata=data[chr==ch]
    ret = MethyLasso:::signal_detection_singlechr(mdata, optimize_lambda2, lambda2, max_lambda2, optimize_hyper, hyper,
                                                 nouter, bf_per_kb, tol_val, max_iter)
    ret$chr=ch
    ret$estimates$chr=ch
    ret
  }
  missing = setdiff(chrs,ret$chr)
  if (length(missing)>0) cat("Warning: the following chromosomes did not work:", missing)
  return(ret)
}

#' Methylation detection using fused lasso, for 1 dataset or several replicates
#'
#' Fast model: normal distribution with identity link and distance-independent modelling
#' 
#' @inheritParams signal_detection
#' 
#' @export

signal_detection_fast = function(data, lambda2=25, tol_val=0.01, max_iter=1, ncores=1,verbose = TRUE) {
  chrs=data[,unique(chr)]
  registerDoParallel(cores=ncores)
  ret = foreach (ch=chrs,.combine=MethyLasso:::combine_ret, .errorhandling = "remove") %dopar% {
    mdata=data[chr==ch]
    if (verbose) cat("Processing chromosome:", ch, "\n")
    ret = MethyLasso:::signal_detection_fast_singlechr(mdata, lambda2, tol_val, max_iter)
    ret$chr=ch
    ret$estimates$chr=ch
    ret
  }
  missing = setdiff(chrs,ret$chr)
  if (length(missing)>0 && verbose) cat("Warning: the following chromosomes did not work:", missing)
  return(ret)
}


#' Methylation difference detection using fused lasso, for 2 datasets (each can
#' have replicates)
#' 
#' beta binomial distribution with logit link
#' 
#' @param ref the ID of the reference condition
#' @inheritParams signal_detection
#' 
#' @export
#' @name difference_detection
difference_detection = function(data, ref, optimize_lambda2=F, lambda2=25, max_lambda2=1000,
                                optimize_hyper=F, hyper=100, nouter=20, bf_per_kb=50, tol_val=0.01,
                                max_iter=30, ncores=1) {
  chrs=data[,unique(chr)]
  registerDoParallel(cores=ncores)
  ret = foreach (ch=chrs,.combine=MethyLasso:::combine_ret, .errorhandling = "remove") %dopar% {
    mdata=data[chr==ch]
    ret = MethyLasso:::difference_detection_singlechr(mdata, ref, optimize_lambda2, lambda2, max_lambda2, optimize_hyper, hyper,
                                                     nouter, bf_per_kb, tol_val, max_iter)
    ret$chr=ch
    ret$estimates$chr=ch
    ret
  }
  missing = setdiff(chrs,ret$chr)
  if (length(missing)>0) cat("Warning: the following chromosomes did not work:", missing)
  return(ret)
}


#' Methylation difference detection using fused lasso, for 2 datasets (each can
#' have replicates)
#' 
#' Fast model: normal distribution with identity link and distance-independent modelling
#' 
#'
#' @param ref the ID of the reference condition
#' @inheritParams signal_detection
#' 
#' @export&& verbose
#' 
difference_detection_fast = function(data, ref, lambda2=25, tol_val=0.01, max_iter=1, ncores=1,verbose = TRUE) {
  chrs=data[,unique(chr)]
  registerDoParallel(cores=ncores)
  ret = foreach (ch=chrs,.combine=MethyLasso:::combine_ret, .errorhandling = "remove") %dopar% {
    mdata=data[chr==ch]
    if (verbose) cat("Processing chromosome:", ch, "\n")
    ret = MethyLasso:::difference_detection_fast_singlechr(mdata, ref, lambda2, tol_val, max_iter)
    ret$chr=ch
    ret$estimates$chr=ch
    ret
  }
  missing = setdiff(chrs,ret$chr)
  if (length(missing)>0 && verbose) cat("Warning: the following chromosomes did not work:", missing)
  return(ret)
}


