#!/opt/R/current/Rscript --vanilla
message("\n\n\n*      *   * * * *   * * * *   *     *   *     *   *            *      * * * *   * * * *   * * * *")
message("* *  * *   *            *      *     *    *   *    *           * *     *         *         *     *")
message("*  **  *   * * *        *      * * * *     * *     *          *   *    * * * *   * * * *   *     *")
message("*      *   *            *      *     *      *      *         * * * *         *         *   *     *")
message("*      *   * * * *      *      *     *     *       * * * *  *       *  * * * *   * * * *   * * * *\n")


## Load libraries
suppressPackageStartupMessages(library(MethyLasso))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(R.utils))

## Retrieve arguments
args=commandArgs(TRUE)


## Help
help <- function(){
    cat("MethyLasso identifies low-methylated regions (LMRs), unmethylated regions (UMRs), DNA methylation valleys (DMVs) and partially methylated domains (PMDs) in a single condition as well as differentially methylated regions (DMRs) between two conditions\n")
    cat("\nUSAGE for two conditions:\nRscript methyLasso.R --n1 [name1] --c1 [cond1_rep1,cond1_rep2,..] --n2 [name2] --c2 [cond2_rep1,cond2_rep2,..] [options]\n")
    cat("\nUSAGE for one condition:\nRscript methyLasso.R --n1 [name1] --c1 [cond1_rep1,cond1_rep2,..] [options]\n")
    cat("\nARGUMENTS:\n")
    cat(" --n1\tName of condition 1 (used in output file names)\n")
    cat(" --c1\tInput file(s) for condition 1 (separated by comma for replicates)\n")
    cat(" --n2\tName of condition 2 (used in output file names and corresponding to reference condition for DMR identification) \n")
    cat(" --c2\tInput file(s) for condition 2 (separated by comma for replicates)\n")
    cat("\nINPUT FORMAT:\n")
    cat("Default: Bismark output file (.gz) with for each CpG: chr / start / end / count of methylated Cs / count of unmethylated Cs\n")
    cat("Alternative: Tab delimited file with for each CpG: chr / start / end, and specify which column contains: count of methylated and unmethylated Cs OR coverage and methylation percentage.\n")
    cat("\n --mC\tColumn number containing count of methylated Cs\n")
    cat(" --uC\tColumn number containing count of unmethylated Cs\n")
    cat("or\n")
    cat(" --cov\tColumn number containing coverage\n")
    cat(" --meth\tColumn number containing percentage of methylation [0-1]\n")
    cat("\nOPTIONS:\n")
    cat(" -c\tMinimum coverage per CpG (default 5)\n")
    cat(" -n\tMinimum number of CpGs in LMRs, UMRS, DMVs and DMRs (default 4)\n")
    cat(" -m\tMaximum mean methylation for PMDs (default 0.7)\n")
    cat(" -d\tMinimum methylation difference for DMRs [0-1] (default 0.1)\n")
    cat(" -p\tP-value threshold for DMRs (default 0.05)\n")
    cat(" -q\tQ-value (FDR) threshold for DMRs (to be used instead of p-value threshold)\n")
    cat(" -r\tReplicate coverage score for CpGs in DMRs (default 0.7)\n")
    cat(" -s\tSkip LMR, UMR, DMV and PMD identification (default FALSE)\n")
    cat(" -t\tNumber of threads to use (default 1)\n")
    cat(" -o\tOutput directory (default current directory)\n") 
    cat(" -f\tOutput figures in pdf format (default TRUE)\n")  
    cat("\nOTHERS OPTIONS:\n")
    cat(" --quiet\tDo not print processing information (default FALSE)\n")
    cat(" --version \tPrint version\n")
    cat(" --help \tPrint help\n")
    cat("\n")
    q()
}


## Default arguments
c = 5
n = 4
m = 0.7
d = 0.1
p = 0.05
q = 1
r = 0.7
t = 1
mC = 4
uC = 5
s = FALSE
f = TRUE
quiet = FALSE


## Save values of each argument
	if(length(args)==0 || !is.na(charmatch("--help",args))|| !is.na(charmatch("-h",args))){
    	help()
	} else {
    	for(ii in 1:length(args)){
		if(grepl("^-",args[ii]) && args[ii] != "-"){
	    	if(ii+1<=length(args) && (!grepl("^-",args[ii+1]) || args[ii+1]=="-")){
	        	assign(gsub("-","",args[ii]),args[ii+1])
	    	} else { assign(gsub("-","",args[ii]),1) }
		}
	}}


s=as.logical(s)
f=as.logical(f)
quiet=as.logical(quiet)



# Package version 
if (!is.na(charmatch("--version", args)) || !is.na(charmatch("-v", args))) {
  package_info <- packageDescription("MethyLasso")
  package_version <- package_info$Version
  cat("MethyLasso version:", package_version, "\n")
  q()
}


# create methylasso output directory 
if (!exists("o")) {
  o = "./"
} else {
  if (!dir.exists(o)) {
    dir.create(o)
  }
}



# Check if we have 1 or 2 conditions
if (exists("c2")) {
  condi <- c(c1, c2)
  nam <- c(n1, n2)
} else {
  condi <- c1
  nam <- n1
}

  
process_replicate <- function(file_path, name, replicate_name) {
  df <- fread(file_path)
  if (exists("cov") && exists("meth")) {
    coverage <- df[[paste0("V", cov)]]
    if (max(df[[paste0("V", meth)]]) > 1) {
      Nmeth <- round(coverage * df[[paste0("V", meth)]] / 100)
      beta <- df[[paste0("V", meth)]] / 100
    } else {
      Nmeth <- round(coverage * df[[paste0("V", meth)]])
      beta <- df[[paste0("V", meth)]]
    }
  } else {
    coverage <- round(df[[paste0("V", mC)]] + df[[paste0("V", uC)]])
    Nmeth <- df[[paste0("V", mC)]]
    beta <- Nmeth / coverage
  }
  
  

  new_df <- data.table(
    dataset = paste0(name, "_", replicate_name),
    condition = name,
    replicate = replicate_name,
    chr = df[["V1"]],
    pos = df[["V2"]],
    coverage = coverage,
    Nmeth = Nmeth,
    beta = beta
  )
  
  new_df[, `:=`(dataset = factor(dataset, unique(dataset)),
                condition = factor(condition, unique(condition)),
                replicate = factor(replicate, unique(replicate)))]

  return(new_df)
}



# Process each condition and replicate using lapply
data_list <- lapply(seq_along(condi), function(i) {
  cond <- condi[i]
  name <- nam[i]

  if (str_detect(cond, ",") == TRUE) {
    files <- str_split_fixed(cond, ",", str_count(cond, ",") + 1)
  } else {
    files <- cond
  }

  if (isFALSE(quiet)) {
    message(paste("Reading condition", name, "with", length(files), "replicate(s) in memory.."))
  }

  replicate_list <- lapply(seq_along(files), function(j) {
    replicate_name <- j
    file_path <- files[j]
    new_df <- process_replicate(file_path, name, replicate_name)
    return(new_df)
  })

  replicate_data <- rbindlist(replicate_list)
  return(replicate_data)
})



# Combine all the data into a single data table
data <- rbindlist(data_list)

# Define a function to generate the plots for a replicate
generate_plots <- function(data, name, replicate_name, plot=f) {
  # Compute mean depth and total CpG count
  meandepth <- mean(data$coverage)
  totcpg <- nrow(data)
  
  # Compute CpG coverage
  cpgcov <- nrow(data[data$coverage >= c,]) * 100 / totcpg
  # Subset the data to select CpGs with coverage >= depth
  a <- which(data$coverage >= c)

 # Set the file name for the PDF
  if (isTRUE(f)) {
  pdf(paste0(o,"/", name, "_", replicate_name, "_plots.pdf"))
  
    # Generate the histograms
    par(mfrow = c(2, 1), bg = "white")
    hist(
      data$coverage[data$coverage < 200], xlim = c(0, 100), breaks = 100,
      main = sprintf("Sequencing Depth Mean=%s x", round(meandepth, digits = 2)),
      xlab = "Read Counts"
    )
    hist(
      data$beta[a] * 100, breaks = 20,
      main = paste0(
        "Distribution of CpG methylation ", c, "X (", round(cpgcov, digits = 2),
        "%)\n Total CpG count = ", totcpg
      ),
      xlab = "DNA methylation (%)"
    )
    null = dev.off()
  }
}

# Generate the plots for each replicate of each condition using lapply
invisible(lapply(data_list, function(x) {
  name <- unique(x$condition)
  replicate_names <- unique(x$replicate)
  
  lapply(replicate_names, function(y) {
    replicate_data <- x[x$replicate == y, ]
    generate_plots(replicate_data, name, y)
  })
}))


## Optionnally filtered by coverage
mindepth = if (!is.null(c)) {
  if (c >= 5) {
    c
  } else {
    c
  }
} else mindepth

data = data[coverage >= mindepth]


 if (isFALSE(s) && exists("c1") && exists("c2")) {
  ### Part 1 and Part2
  # FIT
  if (isFALSE(quiet)) {
    message("\nPart 1: Beginning analysis of the levels of DNA methylation in a single condition\n")
    message("Step 1: Fit of methylation signal\n")
  }

  ret <- MethyLasso:::signal_detection_fast(data, lambda2 = 25, ncores = t, verbose = !quiet)

  # SEGMENT
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of PMDs, LMRs, UMRs and DMVs\n")
  }

  segments <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = m, min_num_cpgs = n)
  seg <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = 1, pmd_std_threshold=0, min_num_cpgs = n)
  # Save file
  for (i in seq_along(data_list)) {
    name <- nam[i]
    # PMD
    pmd <- segments$pmd[condition == name]
    write.table(
      pmd[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = beta, std = std, category = category)],
      file = paste(o, "/", name, "_pmd.tsv", sep = ""),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )

    seg1 <- seg$pmd[condition == name]
    print(seg1)
    # Set the file name for the PDF
    if (isTRUE(f)) {
    	pdf(paste0(o,"/", name, "_pmds_plots.pdf"))
    	par(mfrow = c(1, 1), bg = "white")
    smoothScatter(seg1$beta,seg1$std, ylim=c(0,0.5), xlim=c(0,1), ylab="Standard deviation", xlab="mean DNA methylation", main=paste("MethyLasso segments \n", name))
    abline(v = 0.7, col = rgb(1, 0, 0, alpha = 0.5)) 
    null = dev.off()
    }
	  
    # LMR, UMR and DMV
    lmr_umr_dmv <- segments$lmr_umr_valley[condition == name]
    write.table(
      lmr_umr_dmv[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = beta, std = std, category = category)],
      file = paste(o, "/", name, "_lmr_umr_dmv.tsv", sep = ""),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
  }
	 
ref = n2

  if (isFALSE(quiet)) {
    message("Part 2: Beginning identification of differentially methylated regions (DMRs) across two conditions \n")
    message("Comparison ", n1, " to ", n2, " (ref)")
  } else {
    stop("You need to specify a second condition, name, and reference to find DMRs.")
  }

  # FIT
  if (isFALSE(quiet)) {
    message("\nStep 1: Fit of methylation signal based on both conditions\n")
  }

  ret <- MethyLasso:::difference_detection_fast(data, ref = ref, lambda2 = 25, ncores = t, verbose = !quiet)

  # SEGMENT AND CALL DMRS
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of differentially methylated regions (DMRs)\n")
  }

  diff_call = MethyLasso:::call_differences(data, ret, min.diff = d, pval.cutoff = p, fdr.cutoff = q, cov.score = r * 100, ncores = t, verbose = !quiet, min.num.cpgs = n)

  differences = MethyLasso:::annotate_differences(diff_call, segments)

  write.table(
    differences[, .(chr = chr, start = start, end = end, num.cpgs1 = num.cpgs, num.cpgs2 = num.cpgs.ref, cov.score = coverage.score, meth1 = beta, meth2 = beta.ref, diff = diff, pvalue = pval, FDR = fdr, annotation = switch)],
    file = paste(o, "/", n1, "_vs_", n2, "_dmrs.tsv", sep = ""),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )

  # Create DMR plot
  if (isTRUE(f)) {
    # DMR scatterplot
    pdf(paste0(o, "/", n1, "_vs_", n2, "_dmrs_plot.pdf", sep = ""))
    par(mfrow = c(1, 1), bg = "white")
    smoothScatter(
      diff_call$beta * 100, diff_call$beta.ref * 100,
      xlab = paste0(n1, " DNA methylation (%)"),
      ylab = paste0(n2, " DNA methylation (%)"),
      main = paste0("DNA methylation ", n1, " vs ", n2)
    )
    null = dev.off()
  }
} else if (isTRUE(s) && exists("c1") && exists("c2")) {
  ### Only Part2
  # Create name n which is not the reference name
  ref = n2

  if (isFALSE(quiet)) {
    message("\nPart 2: Beginning identification of differentially methylated regions (DMRs) across two conditions \n")
    message("Comparison ", n1, " to ", n2, " (ref)")
  } else {
    stop("You need to specify a second condition, name, and reference to find DMRs.")
  }

  # FIT
  if (isFALSE(quiet)) {
    message("\nStep 1: Fit of methylation signal based on both conditions\n")
  }

  ret <- MethyLasso:::difference_detection_fast(data, ref = ref, lambda2 = 25, ncores = t, verbose = !quiet)

# SEGMENT AND CALL DMRS
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of differentially methylated regions (DMRs)\n")
  }

  diff_call = MethyLasso:::call_differences(data, ret, min.diff = d, pval.cutoff = p, fdr.cutoff = q, min.num.cpgs = n, cov.score = r * 100, ncores = t, verbose = !quiet)

  write.table(
    diff_call[, .(chr = chr, start = start, end = end, num.cpgs1 = num.cpgs, num.cpgs2 = num.cpgs.ref, cov.score = coverage.score, meth1 = beta, meth2 = beta.ref, diff = diff, pvalue = pval, FDR = fdr)],
    file = paste(o, "/", n1, "_vs_", n2, "_dmrs.tsv", sep = ""),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
  # Create DMR plot
  if (isTRUE(f)) {
    # DMR scatterplot
    pdf(paste0(o, "/", n1, "_vs_", n2, "_dmrs_plot.pdf", sep = ""))
    par(mfrow = c(1, 1), bg = "white")
    smoothScatter(
      diff_call$beta * 100, diff_call$beta.ref * 100,
      xlab = paste0(n1, " DNA methylation (%)"),
      ylab = paste0(n2, " DNA methylation (%)"),
      main = paste0("DNA methylation ", n1, " vs ", n2)
    )
    null = dev.off()
  }

} else if ( exists("c1") && !exists("c2")|| exists("c2")  && !exists("c1")) {
  ### Only Part1
  if (isFALSE(quiet)) {
    message("\nPart 1: Beginning analysis of the levels of DNA methylation in a single condition\n")
    message("Step 1: Fit of methylation signal\n")
  }
  # FIT
  ret <- MethyLasso:::signal_detection_fast(data, lambda2 = 25, ncores = t, verbose = !quiet)

  # SEGMENT
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of PMDs, LMRs, UMRs and DMVs\n")
  }

segments <- MethyLasso:::segment_methylation(data, ret, ncores = t, min_num_cpgs = n, pmd_max_beta = m)
seg <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = 1, pmd_std_threshold=0, min_num_cpgs = n)
	 
  # Save file
  for (i in seq_along(data_list)) {
    name <- nam[i]
    # PMD
    pmd <- segments$pmd[condition == name]
    write.table(
      pmd[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = beta, std = std, category = category)],
      file = paste(o, "/", name, "_pmd.tsv", sep = ""),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
	  
    seg1 <- seg$pmd[condition == name]
    print(seg1)
    # Set the file name for the PDF
    if (isTRUE(f)) {
    	pdf(paste0(o,"/", name, "_pmds_plots.pdf"))
    	par(mfrow = c(1, 1), bg = "white")
    smoothScatter(seg1$beta,seg1$std, ylim=c(0,0.5), xlim=c(0,1), ylab="Standard deviation", xlab="mean DNA methylation", main=paste("MethyLasso segments \n ",name))
    abline(v = 0.7, col = rgb(1, 0, 0, alpha = 0.5)) 
    null = dev.off()
    }
	
    # LMR, UMR and DMV
    lmr_umr_dmv <- segments$lmr_umr_valley[condition == name]
    write.table(
      lmr_umr_dmv[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = beta, std = std, category = category)],
      file = paste(o, "/", name, "_lmr_umr_dmv.tsv", sep = ""),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
  }

} else {
  stop("See help ")
}
