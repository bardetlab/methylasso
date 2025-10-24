#!/opt/R/current/Rscript --vanilla
message("\n\n\n*      *   * * * *   * * * *   *     *   *     *   *            *      * * * *   * * * *   * * * *")
message("* *  * *   *            *      *     *    *   *    *           * *     *         *         *     *")
message("*  **  *   * * *        *      * * * *     * *     *          *   *    * * * *   * * * *   *     *")
message("*      *   *            *      *     *      *      *         * * * *         *         *   *     *")
message("*      *   * * * *      *      *     *     *       * * * *  *       *  * * * *   * * * *   * * * *\n")

options(scipen = 999)

## Load libraries
suppressPackageStartupMessages({
  library(MethyLasso)
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(stringr)
})

# Load R.utils separately to suppress its specific warning
suppressWarnings(suppressPackageStartupMessages(library(R.utils)))

## Function to display memory usage
show_memory_usage <- function(step_name) {
  gc_result <- invisible(gc())  # Suppress the Ncells/Vcells output
  # gc() returns values in Mb in columns 2 and 4
  used_gb <- sum(gc_result[, 2]) / 1024  # Column 2 is "used (Mb)"
  max_used_gb <- sum(gc_result[, 4]) / 1024  # Column 4 is "max used (Mb)"
  message(sprintf("\n[Memory Usage - %s] Current: %.1f GB | Peak: %.1f GB\n", 
                  step_name, used_gb, max_used_gb))
}

## Retrieve arguments
args <- commandArgs(TRUE)

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
  cat("\nADDITIONAL OPTIONS FOR ADVANCED CONFIGURATIONS:\n")
  cat(" --umrlmr_lambda2\tSegmentation value for UMR/LMR identification (default 25)\n")
  cat(" --tol_val\t\tTolerance value for optimization (default 0.01)\n")
  cat(" --pmd_lambda2\t\tSegmentation value for PMD identification (default 1000)\n")
  cat(" --pmd_std_threshold\tStandard deviation threshold for PMD (default 0.15)\n")
  cat(" --umr_std_threshold\tStandard deviation threshold for UMR (default 0.1)\n")
  cat(" --umr_max_beta\t\tMaximum methylation for UMR (default 0.1)\n")
  cat(" --lmr_max_beta\t\tMaximum methylation for LMR (default 0.5)\n")
  cat(" --valley_max_beta\tMaximum methylation for valleys (default 0.1)\n")
  cat(" --pmd_valley_min_width\tMinimum width for PMD valleys (default 5000)\n")
  cat(" --max_distance\t\tMaximum distance for merging segments (default 500)\n")
  cat(" --min_width\t\tMinimum width for segments (default 30)\n")
  cat(" --umr_large_width\tLarge width threshold for UMRs (default 500)\n")
  cat(" --umr_large_density\tDensity threshold for large UMRs (default 0.03)\n")
  cat(" --dmr_lambda2\t\tSegmentation value for DMR identification (default 25)\n")
  cat(" --min.overlap.pc\tMinimum overlap percentage for DMRs (default 10)\n")
  cat(" --flank.width\t\tMinimum width of flanking segments (default 300)\n")
  cat(" --flank.dist.max\tMaximum distance of UMR/LMR to flanking segment (default 5000)\n")
  cat(" --flank.num.cpgs\tConsider as small UMR/LMRs those which have num.cpgs smaller than this threshold (default 10)\n")
  cat(" --flank.beta.min\tFor small UMR/LMRs, require that flanking segments have mean methylation above this threshold (default 0.5)\n")
  cat(" --split.pmds\t\tWhether to split PMDs which contain UMRs (default TRUE)\n")
  cat(" --merge.pmds\t\tWhether to merge adjacent PMDs into a single region (default TRUE)\n")
  cat(" --drop\t\t\tWhether to drop verbose columns and unannotated segments (default TRUE)\n")
  cat("\n")
  cat("\nOTHERS OPTIONS:\n")
  cat(" --quiet\tDo not print processing information (default FALSE)\n")
  cat(" --version \tPrint version\n")
  cat(" --help \tPrint help\n")
  cat("\n")
  q()
}

## Default arguments
mC = 5 # Column number containing count of methylated Cs
uC = 6 # Column number containing count of unmethylated Cs
c = 5 # Minimum coverage per CpG
n = 4 # Minimum number of CpGs in LMRs, UMRS, DMVs and DMRs
m = 0.7 # Maximum mean methylation for PMDs
d = 0.1 # Minimum methylation difference for DMRs
p = 0.05 # P-value threshold for DMRs
q = 1 # Q-value (FDR) threshold for DMRs
r = 0.7 # Replicate coverage score for CpGs in DMRs
s = FALSE # Skip LMR, UMR, DMV and PMD identification 
t = 1 # Number of threads to use
f = TRUE # Output figures in pdf format
quiet = FALSE # Do not print processing information
umrlmr_lambda2 = 25 # Segmentation value for UMR/LMR identification
tol_val = 0.01 # Tolerance value for optimization
pmd_lambda2 = 1000 # Segmentation value for PMD identification
pmd_std_threshold = 0.15 # Standard deviation threshold for PMD
umr_std_threshold = 0.1 # Standard deviation threshold for UMR
umr_max_beta = 0.1 # Maximum mean methylation for UMR
lmr_max_beta = 0.5 # Maximum mean methylation for LMR
valley_max_beta = 0.1 # Maximum beta value for valleys
pmd_valley_min_width = 5000 # Minimum width for PMD valleys
max_distance = 500 # Maximum distance for merging segments
min_width = 30 # Minimum width for segments
umr_large_width = 500 # Large width threshold for UMRs
umr_large_density = 0.03 # Density threshold for large UMRs
dmr_lambda2 = 25 #  Segmentation value for DMR identification
min.overlap.pc = 10 # Minimum overlap percentage for DMRs
flank.width = 300 # Minimum width of flanking segments
flank.dist.max = 5000 # Maximum distance of UMR/LMR to flanking segment
flank.num.cpgs = 10 # Threshold for considering small UMR/LMRs
flank.beta.min = 0.5 # Minimum methylation for flanking segments of small UMR/LMRs
split.pmds = TRUE # Whether to split PMDs which contain UMRs
merge.pmds = TRUE # Whether to merge adjacent PMDs
drop = TRUE # Whether to drop verbose columns and unannotated segments

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
  }
}

# Convert to numeric or logical 
mC = as.numeric(mC)
uC = as.numeric(uC)
c = as.numeric(c)
n = as.numeric(n)
m = as.numeric(m)
d = as.numeric(d)
p = as.numeric(p)
q = as.numeric(q)
r = as.numeric(r)
umrlmr_lambda2 = as.numeric(umrlmr_lambda2)
tol_val = as.numeric(tol_val)
pmd_lambda2 = as.numeric(pmd_lambda2)
pmd_std_threshold = as.numeric(pmd_std_threshold)
umr_std_threshold = as.numeric(umr_std_threshold)
umr_max_beta = as.numeric(umr_max_beta )
lmr_max_beta = as.numeric(lmr_max_beta)
valley_max_beta = as.numeric(valley_max_beta)
pmd_valley_min_width = as.numeric(pmd_valley_min_width)
max_distance = as.numeric(max_distance)
min_width = as.numeric(min_width)
umr_large_width = as.numeric(umr_large_width)
umr_large_density = as.numeric(umr_large_density)
dmr_lambda2 = as.numeric(dmr_lambda2)
min.overlap.pc = as.numeric(min.overlap.pc)
flank.width = as.numeric(flank.width)
flank.dist.max = as.numeric(flank.dist.max)
flank.num.cpgs = as.numeric(flank.num.cpgs)
flank.beta.min = as.numeric(flank.beta.min)

s = as.logical(s)
f = as.logical(f)
quiet = as.logical(quiet)
split.pmds = as.logical(split.pmds)
merge.pmds = as.logical(merge.pmds)
drop = as.logical(drop)

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


# Process data in chunks to reduce memory usage
process_condition <- function(cond, name) {
  if (str_detect(cond, ",")) {
    files <- str_split_fixed(cond, ",", str_count(cond, ",") + 1)
  } else {
    files <- cond
  }
  
  if (!quiet) {
    message(paste("Processing condition", name, "with", length(files), "replicate(s)..."))
  }
  
  # Process replicates sequentially to save memory
  result_list <- vector("list", length(files))
  for (j in seq_along(files)) {
    result_list[[j]] <- process_replicate(files[j], name, j)
    # Generate plots immediately and free memory
    if (f) {
      generate_plots(result_list[[j]], name, j)
    }
    invisible(gc())  # Force garbage collection
  }
  
  # Combine results
  replicate_data <- rbindlist(result_list)
  rm(result_list)
  invisible(gc())
  
  return(replicate_data)
}

# Process each file with memory efficiency
process_replicate <- function(file_path, name, replicate_name) {
  # Read the first row to check if there's a header
  first_row <- fread(file_path, nrows = 1, header = FALSE)
  
  # Determine if the file has a header by checking if the first row contains character data
  has_header <- any(sapply(first_row, is.character))
  
  # Read only necessary columns by index
  if (exists("cov") && exists("meth")) {
    # Use column indices instead of names
    df <- fread(file_path, select = c(1, 2, as.integer(cov), as.integer(meth)), 
                header = has_header, verbose = FALSE, showProgress = !quiet)
  } else {
    # Use column indices instead of names
    df <- fread(file_path, select = c(1, 2, as.integer(mC), as.integer(uC)), 
                header = has_header, verbose = FALSE, showProgress = !quiet)
  }
  
  # Process data
  if (exists("cov") && exists("meth")) {
    coverage <- df[[3]]
    if (max(df[[4]]) > 1) {
      Nmeth <- round(coverage * df[[4]] / 100)
      beta <- df[[4]] / 100
    } else {
      Nmeth <- round(coverage * df[[4]])
      beta <- df[[4]]
    }
  } else {
    coverage <- round(df[[3]] + df[[4]])
    Nmeth <- df[[3]]
    beta <- Nmeth / coverage
  }
  
  # Create efficient data.table
  new_df <- data.table(
    dataset = paste0(name, "_", replicate_name),
    condition = name,
    replicate = replicate_name,
    chr = df[[1]],
    pos = df[[2]],
    coverage = coverage,
    Nmeth = Nmeth,
    beta = beta
  )
  
  # Clean up
  rm(df, coverage, Nmeth, beta)
  invisible(gc())
  
  return(new_df)
}

# Generate plots for a replicate
generate_plots <- function(data, name, replicate_name, plot = f) {
  if (!plot) return()
  
  # Compute mean depth and total CpG count
  meandepth <- mean(data$coverage)
  totcpg <- nrow(data)
  
  # Compute CpG coverage
  cpgcov <- nrow(data[data$coverage >= c, ]) * 100 / totcpg
  # Subset the data to select CpGs with coverage >= depth
  a <- which(data$coverage >= c)
  
  # Set the file name for the PDF
  pdf(paste0(o, "/", name, "_", replicate_name, "_plots.pdf"))
  
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
  invisible(dev.off())  # Use invisible() to suppress the null device message
}

# Process all conditions with memory efficiency
data_list <- vector("list", length(condi))
for (i in seq_along(condi)) {
  data_list[[i]] <- process_condition(condi[i], nam[i])
  invisible(gc())
}

# Combine data and set factors
data <- rbindlist(data_list)
rm(data_list)
invisible(gc())


data[, `:=`(
  dataset = factor(dataset, unique(dataset)),
  condition = factor(condition, unique(condition)),
  replicate = factor(replicate, unique(replicate))
)]

# Filter by coverage
mindepth = if (!is.null(c)) {
  if (c >= 5) {
    c
  } else {
    c
  }
} else mindepth

data <- data[coverage >= mindepth]
invisible(gc())

# Rest of the processing remains the same but with additional invisible(gc()) calls
if (isFALSE(s) && exists("c1") && exists("c2")) {
  ### Part 1 and Part2
  # FIT
  if (isFALSE(quiet)) {
    message("\nPart 1: Beginning analysis of the levels of DNA methylation in a single condition\n")
    message("Step 1: Fit of methylation signal\n")
  }
  
  ret <- MethyLasso:::signal_detection_fast(data, lambda2 = umrlmr_lambda2, tol_val = tol_val, ncores = t, verbose = !quiet)
  invisible(gc())
  
  # SEGMENT
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of PMDs, LMRs, UMRs and DMVs\n")
  }
  
  segments <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = m, min_num_cpgs = n, 
                                               pmd_lambda2 = pmd_lambda2, pmd_std_threshold = pmd_std_threshold, 
                                               umr_std_threshold = umr_std_threshold, umr_max_beta = umr_max_beta, 
                                               lmr_max_beta = lmr_max_beta, valley_max_beta = valley_max_beta, 
                                               pmd_valley_min_width = pmd_valley_min_width, max_distance = max_distance, 
                                               min_width = min_width, umr_large_width = umr_large_width, 
                                               umr_large_density = umr_large_density, flank.width = flank.width, 
                                               flank.dist.max = flank.dist.max, flank.num.cpgs = flank.num.cpgs, 
                                               flank.beta.min = flank.beta.min, split.pmds = split.pmds,
                                               merge.pmds = merge.pmds, drop = drop, tol_val = tol_val)
  invisible(gc())
  
  #to obtain all segments for plot
  seg <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = 1, pmd_std_threshold = 0, 
                                          min_num_cpgs = n, pmd_lambda2 = pmd_lambda2, umr_std_threshold = umr_std_threshold,
                                          umr_max_beta = umr_max_beta, lmr_max_beta = lmr_max_beta, 
                                          valley_max_beta = valley_max_beta, pmd_valley_min_width = pmd_valley_min_width, 
                                          max_distance = max_distance, min_width = min_width, umr_large_width = umr_large_width, 
                                          umr_large_density = umr_large_density, flank.width = flank.width, 
                                          flank.dist.max = flank.dist.max, flank.num.cpgs = flank.num.cpgs, 
                                          flank.beta.min = flank.beta.min, split.pmds = split.pmds,
                                          merge.pmds = merge.pmds, drop = drop, tol_val = tol_val)
  invisible(gc())
  
  # Save file
  for (i in seq_along(nam)) {
    name <- nam[i]
    # PMD
    pmd <- segments$pmd[condition == name]
    write.table(
      pmd[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = round(beta*100, 2), std = round(std,3), category = category)],
      file = paste0(o, "/", name, "_pmd.tsv"),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
    
    seg1 <- seg$pmd[condition == name]
    
    # Set the file name for the PDF
    if (f) {
      pdf(paste0(o, "/", name, "_pmds_plots.pdf"))
      par(mfrow = c(1, 1), bg = "white")
      smoothScatter(seg1$beta*100, seg1$std, ylim = c(0, 0.5), xlim = c(0, 100), 
                    ylab = "Standard deviation", xlab = "DNA methylation (%)", 
                    main = paste("MethyLasso segments \n", name))
      abline(v = 70, col = rgb(1, 0, 0, alpha = 0.5)) 
      invisible(dev.off())  # Use invisible() to suppress the null device message
    }
    
    # LMR, UMR and DMV
    lmr_umr_dmv <- segments$lmr_umr_valley[condition == name]
    write.table(
      lmr_umr_dmv[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = round(beta*100, 2), std = round(std,3), category = category)],
      file = paste0(o, "/", name, "_lmr_umr_dmv.tsv"),
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
  
  show_memory_usage("Before Step 1: Fit of methylation signal")
  ret <- MethyLasso:::difference_detection_fast(data, ref = ref, lambda2 = dmr_lambda2, ncores = t, verbose = !quiet)
  invisible(gc())
  show_memory_usage("After Step 1: Fit of methylation signal")
  
  # SEGMENT AND CALL DMRS
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of differentially methylated regions (DMRs)\n")
  }
  
  show_memory_usage("Before Step 2: DMR identification")
  diff_call = MethyLasso:::call_differences(data, ret, min.diff = d, pval.cutoff = p, fdr.cutoff = q, 
                                            cov.score = r * 100, ncores = t, verbose = !quiet, 
                                            min.num.cpgs = n, tol_val = tol_val)
  invisible(gc())
  show_memory_usage("After Step 2: DMR identification")
  
  differences = MethyLasso:::annotate_differences(diff_call, segments, min.overlap.pc = min.overlap.pc)
  invisible(gc())
  
  write.table(
    differences[, .(chr = chr, start = start, end = end, num.cpgs1 = num.cpgs, num.cpgs2 = num.cpgs.ref, 
                    cov.score = round(coverage.score, 3), meth1 = round(beta * 100, 2), 
                    meth2 = round(beta.ref * 100, 2), diff = round(diff * 100, 2), 
                    pvalue = pval, FDR = fdr, annotation = switch)],
    file = paste0(o, "/", n1, "_vs_", n2, "_dmrs.tsv"),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
  
  # Create DMR plot
  if (f) {
    # DMR scatterplot
    pdf(paste0(o, "/", n1, "_vs_", n2, "_dmrs_plot.pdf"))
    par(mfrow = c(1, 1), bg = "white")
    smoothScatter(
      diff_call$beta * 100, diff_call$beta.ref * 100,
      xlab = paste0(n1, " DNA methylation (%)"),
      ylab = paste0(n2, " DNA methylation (%)"),
      main = paste0("DNA methylation ", n1, " vs ", n2)
    )
    invisible(dev.off())  # Use invisible() to suppress the null device message
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
  
  show_memory_usage("Before Step 1: Fit of methylation signal")
  ret <- MethyLasso:::difference_detection_fast(data, ref = ref, lambda2 = dmr_lambda2, ncores = t, verbose = !quiet)
  invisible(gc())
  show_memory_usage("After Step 1: Fit of methylation signal")
  
  # SEGMENT AND CALL DMRS
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of differentially methylated regions (DMRs)\n")
  }
  
  show_memory_usage("Before Step 2: DMR identification")
  diff_call = MethyLasso:::call_differences(data, ret, min.diff = d, pval.cutoff = p, fdr.cutoff = q, 
                                            cov.score = r * 100, ncores = t, verbose = !quiet, 
                                            min.num.cpgs = n, tol_val = tol_val)
  invisible(gc())
  show_memory_usage("After Step 2: DMR identification")
  
  write.table(
    diff_call[, .(chr = chr, start = start, end = end, num.cpgs1 = num.cpgs, num.cpgs2 = num.cpgs.ref, 
                  cov.score = round(coverage.score, 3), meth1 = round(beta * 100, 2), 
                  meth2 = round(beta.ref * 100, 2), diff = round(diff * 100, 2), 
                  pvalue = round(pval, 4), FDR = round(fdr, 4))],
    file = paste0(o, "/", n1, "_vs_", n2, "_dmrs.tsv"),
    quote = FALSE, row.names = FALSE, sep = "\t"
  )
  
  # Create DMR plot
  if (f) {
    # DMR scatterplot
    pdf(paste0(o, "/", n1, "_vs_", n2, "_dmrs_plot.pdf"))
    par(mfrow = c(1, 1), bg = "white")
    smoothScatter(
      diff_call$beta * 100, diff_call$beta.ref * 100,
      xlab = paste0(n1, " DNA methylation (%)"),
      ylab = paste0(n2, " DNA methylation (%)"),
      main = paste0("DNA methylation ", n1, " vs ", n2)
    )
    invisible(dev.off())  # Use invisible() to suppress the null device message
  }
  
} else if (exists("c1") && !exists("c2") || exists("c2") && !exists("c1")) {
  ### Only part 1
  if (isFALSE(quiet)) {
    message("\nPart 1: Beginning analysis of the levels of DNA methylation in a single condition\n")
    message("Step 1: Fit of methylation signal\n")
  }
  show_memory_usage("Before Step 1: Fit of methylation signal")  # Added memory check
  # FIT
  ret <- MethyLasso:::signal_detection_fast(data, lambda2 = umrlmr_lambda2, tol_val = tol_val, ncores = t, verbose = !quiet)
  invisible(gc())
  show_memory_usage("After Step 1: Fit of methylation signal")  # Added memory check
  
  # SEGMENT
  if (isFALSE(quiet)) {
    message("\nStep 2: Segmentation and identification of PMDs, LMRs, UMRs and DMVs\n")
  }
  show_memory_usage("Before Step 2: Segmentation")  # Added memory check
  segments <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = m, min_num_cpgs = n, 
                                               pmd_lambda2 = pmd_lambda2, pmd_std_threshold = pmd_std_threshold, 
                                               umr_std_threshold = umr_std_threshold, umr_max_beta = umr_max_beta, 
                                               lmr_max_beta = lmr_max_beta, valley_max_beta = valley_max_beta, 
                                               pmd_valley_min_width = pmd_valley_min_width, max_distance = max_distance, 
                                               min_width = min_width, umr_large_width = umr_large_width, 
                                               umr_large_density = umr_large_density, flank.width = flank.width, 
                                               flank.dist.max = flank.dist.max, flank.num.cpgs = flank.num.cpgs, 
                                               flank.beta.min = flank.beta.min, split.pmds = split.pmds,
                                               merge.pmds = merge.pmds, drop = drop, tol_val = tol_val)
  invisible(gc())
  show_memory_usage("After Step 2: Segmentation")  # Added memory check
  
  #to obtain all segments for plot
  seg <- MethyLasso:::segment_methylation(data, ret, ncores = t, pmd_max_beta = 1, pmd_std_threshold = 0, 
                                          min_num_cpgs = n, pmd_lambda2 = pmd_lambda2, umr_std_threshold = umr_std_threshold,
                                          umr_max_beta = umr_max_beta, lmr_max_beta = lmr_max_beta, 
                                          valley_max_beta = valley_max_beta, pmd_valley_min_width = pmd_valley_min_width, 
                                          max_distance = max_distance, min_width = min_width, umr_large_width = umr_large_width, 
                                          umr_large_density = umr_large_density, flank.width = flank.width, 
                                          flank.dist.max = flank.dist.max, flank.num.cpgs = flank.num.cpgs, 
                                          flank.beta.min = flank.beta.min, split.pmds = split.pmds,
                                          merge.pmds = merge.pmds, drop = drop, tol_val = tol_val)
  invisible(gc())
  show_memory_usage("After obtaining all segments for plotting")  # Added memory check
  
  # Save file
  for (i in seq_along(nam)) {
    name <- nam[i]
    # PMD
    pmd <- segments$pmd[condition == name]
    write.table(
      pmd[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = round(beta*100, 2), std = round(std,3), category = category)],
      file = paste0(o, "/", name, "_pmd.tsv"),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
    
    seg1 <- seg$pmd[condition == name]
    
    # Set the file name for the PDF
    if (f) {
      pdf(paste0(o, "/", name, "_pmds_plots.pdf"))
      par(mfrow = c(1, 1), bg = "white")
      smoothScatter(seg1$beta*100, seg1$std, ylim = c(0, 0.5), xlim = c(0, 100), 
                    ylab = "Standard deviation", xlab = "DNA methylation (%)", 
                    main = paste("MethyLasso segments \n ", name))
      abline(v = 70, col = rgb(1, 0, 0, alpha = 0.5)) 
      invisible(dev.off())
    }
    
    # LMR, UMR and DMV
    lmr_umr_dmv <- segments$lmr_umr_valley[condition == name]
    write.table(
      lmr_umr_dmv[, .(chr = chr, start = start, end = end, num.cpgs = num.cpgs, meth = round(beta*100, 2), std = round(std,3), category = category)],
      file = paste0(o, "/", name, "_lmr_umr_dmv.tsv"),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
  }
  
} else {
  stop("See help ")
}
