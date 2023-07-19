	#' @include MethyLasso.R
NULL

#' Segment methylation data (single chromosome version)
#' 
#' @param data the input data frame
#' @param ret the output of signal_detection
#' @param std_lambda2 the value of the fusion penalty for PMD detection (default is 1000)
#' @param pmd_std_threshold the value of the threshold on the std (default 0.15) used for PMDs
#' @param umr_std_threshold the value of the threshold on the std (default 0.1)  use for UMRs
#' @param umr_max_beta segments with mean methylation below this threshold are UMRs (default 0.1)
#' @param lmr_max_beta segments with mean methylation below this threshold (default 0.5) and above umr_max_beta are LMRs
#' @param pmd_max_beta the maximum methylation value a segment group can have to be considered a PMD (default 0.70)
#' @param valley_max_beta Valleys must have mean methylation below this threshold (default 0.1)
#' @param pmd_valley_min_width Valleys and PMDs must have at least this size (default 0)
#' @param max_distance Within a segment, the distance between successive CpGs must never be larger than this distance,
#'   and its mean must also be below it. (default )
#' @param min_num_cpgs minimum number of CpGs to classify as LMR/UMR (default 4)
#' @param min_width minimum width to classify as LMR/UMR (default 10)
#' @param flank.width minimum width of flanking segments (default 300)
#' @param flank.dist.max maximum distance of UMR/LMR to flanking segment (default 0)
#' @param flank.num.cpgs Consider as small UMR/LMRs those which have num.cpgs smaller than this threshold (default 10) 
#' @param flank.beta.min For small UMR/LMRs, require that flanking segments have mean methylation above this threshold (default 0.5)
#' @param umr_large_width large UMRs are accepted when larger than this (default )
#' @param umr_large_density large UMRs are accepted when more dense than this (default 0.03)
#' @param split.pmds whether to split PMDs which contain UMRs (default TRUE)
#' @param merge.pmds whether to merge adjacent PMDs into a single region (default TRUE)
#' @param drop whether to drop verbose columns and unannotated segments (default TRUE)
#' @param tol_val fused values below tol_val from each other are considered equal (default 0.01)
#' 
segment_methylation_singlechr = function(data, ret, std_lambda2=1000, pmd_std_threshold=0.15, umr_std_threshold=0.1,
                                         umr_max_beta=0.1, lmr_max_beta=0.5, pmd_max_beta=0.75, valley_max_beta=0.1, pmd_valley_min_width=5000,
                                         max_distance=500, min_num_cpgs=4, min_width=30, flank.width=300, flank.dist.max=5000,
                                         flank.num.cpgs=10, flank.beta.min=0.5, umr_large_width=500, umr_large_density=0.03,
                                         split.pmds=T, merge.pmds=T, drop=T, tol_val=0.01) {

  # call segments and bin original data along them
  sig.df = ret$estimates
  mdata = copy(as.data.table(data))
  if (ret$detection.type != "signal") {
    # Segment_methylation called on a difference calculation, will return the segmentation of the reference dataset
    sig.df = sig.df[sig.df$estimate_id==1,]
    mdata = mdata[condition==ret$ref.cond]
  }



  mdata[,coverage:=sum(coverage),by=c("condition","chr","pos")]
  mdata[,Nmeth:=sum(Nmeth),by=c("condition","chr","pos")]
  mdata[,beta:=Nmeth/coverage]
  mdata=unique(subset(mdata,,c(condition,chr,pos,coverage,Nmeth,beta)))

  segmentation = as.data.table(MethyLasso:::segment_methylation_helper(sig.df, tol_val = tol_val))
  stopifnot(segmentation[,uniqueN(estimate_id)] == mdata[,uniqueN(condition)])
  setkey(segmentation,estimate_id,start,end)

  # Collect methylation data per segment
  mdata[,c("estimate_id","start","end"):=list(as.integer(condition),pos,pos+1)]
  setkey(mdata,estimate_id,start,end)
  sg <- segmentation[start <= end]
  ov = foverlaps(mdata,sg)

  # Split segments containing gaps
  ov[,is.oversize:=pos-shift(pos)>max_distance,keyby=c("condition","estimate_id","start","end")]
  ov[is.na(is.oversize),is.oversize:=F]
  ov[,segment_id:=segment_id+cumsum(is.oversize),keyby=c("condition","estimate_id")]

  # Shrink boundaries to span of CpGs
  ov[,c("start","end"):=list(min(pos),max(pos)+2),keyby=c("condition","estimate_id","segment_id")]
  ov[,width:=end-start+1]

  # Compute averages on each segment
  seg_call = ov[,.(beta=mean(Nmeth/coverage),coverage=sum(coverage),pseudoweight=pseudoweight[1],
                   num.cpgs=length(unique(pos)),std=sd(Nmeth/coverage)),
                 keyby=c("condition","estimate_id","segment_id","start","end","width")]

  seg_call[is.na(std),std:=0]
  setkeyv(seg_call,c("condition","start","end"))

  # Annotate LMRs/UMRs
  seg_call[,is.admissible:=width>min_width&width/(num.cpgs-1)<max_distance&num.cpgs>=min_num_cpgs]
  seg_call[is.admissible==T,is.V:=beta<shift(beta)&beta<shift(beta,type = "lead"),by=c("condition")]
  seg_call[is.V==T|(width>flank.width&is.admissible==T), is.flank:=is.V!=T & (shift(is.V)==T | shift(is.V,type="lead")==T)]
  seg_call[is.V==T|is.flank==T,flank.dist:=pmax(shift(start,type="lead")-end,start-shift(end)),by=c("condition")]
  seg_call[is.V==T|is.flank==T,flank.beta:=pmin(shift(beta,type="lead"),shift(beta)),by=c("condition")]
  seg_call[is.V==T|is.flank==T,is.sep:=pmin(abs(beta-shift(beta)),abs(beta-shift(beta,type="lead")))>2*std,by=c("condition")]
  seg_call[is.V==T&is.sep==T&flank.dist<flank.dist.max&beta<lmr_max_beta,category:="LMR"]
  seg_call[category=="LMR"&num.cpgs<flank.num.cpgs&flank.beta<=flank.beta.min,category:=NA]
  seg_call[beta<umr_max_beta&std<umr_std_threshold&category=="LMR",category:="UMR"]
  seg_call[beta<umr_max_beta&std<umr_std_threshold&(is.V==T&width>umr_large_width&num.cpgs/width>umr_large_density),category:="UMR"]
  
  # Fuse std values weighted by width
  seg_call[,fused_std:=MethyLasso:::weighted_1d_flsa(std,width,std_lambda2),by=c("condition")]
  pmd.df = as.data.frame(seg_call[,.(estimate_id,start,beta=fused_std,betahat=std,pseudoweight)])
  std_call = as.data.table(MethyLasso:::segment_methylation_helper(pmd.df,tol_val=tol_val))
  setnames(std_call,c("beta","betahat"),c("fused_std","std"))
  setkey(std_call,estimate_id,start,end)
  split_call=std_call

  # Annotate Valleys
  # find low methylation regions
  seg_call[,is.high:=beta>valley_max_beta]
  seg_call[,follows.large.gap:=start-shift(end)>pmd_valley_min_width,by=condition]
  seg_call[is.na(follows.large.gap),follows.large.gap:=F]
  seg_call[,low.id:=cumsum(is.high)+cumsum(follows.large.gap),by=c("condition")]
  seg_call[is.high==F,c("low.id.width","low.id.beta"):=list(max(end)-min(start)+1,mean(beta)),by=c("low.id","condition")]
  new_valleys=seg_call[is.high==F&low.id.beta<=valley_max_beta][
    ,.(start=min(start),end=max(end),width=low.id.width[1],beta=low.id.beta[1],coverage=sum(coverage),
       pseudoweight=sum(pseudoweight),num.cpgs=sum(num.cpgs),category="DMV"),by=c("low.id","low.id.width","low.id.beta","condition")]
  new_valleys[,category:=ifelse(low.id.width>pmd_valley_min_width,"DMV","low")]
  setkey(new_valleys,condition,start,end)
  # add UMRs which are not within these regions
  ov = foverlaps(seg_call[category=="UMR"],new_valleys,which=T,mult="all")
  new_valleys=rbind(new_valleys,seg_call[category=="UMR",.(low.id,low.id.width,low.id.beta,condition,start,end,width=low.id.width,
                                                           beta=low.id.beta,coverage,pseudoweight,num.cpgs,category)][ov[is.na(yid),xid]])
  setkey(new_valleys,condition,start,end)
  
  # Merge consecutive Valleys
  new_valleys[,is.sep.valley:=start-shift(end)>pmd_valley_min_width,by=condition]
  new_valleys[is.na(is.sep.valley),is.sep.valley:=F]
  new_valleys[,valley.id:=cumsum(is.sep.valley),by=condition]
  new_valleys=new_valleys[
    ,.(start=min(start),end=max(end),width=max(end)-min(start)+1,beta=mean(beta),coverage=sum(coverage),
       pseudoweight=sum(pseudoweight),num.cpgs=sum(num.cpgs),contains.valley="DMV"%in%category),by=c("valley.id","condition")]
  new_valleys[!(width>pmd_valley_min_width&width/(num.cpgs-1)<max_distance&num.cpgs>=min_num_cpgs),contains.valley:=F]
  new_valleys=new_valleys[contains.valley==T] 
  new_valleys[,c("valley.id","contains.valley"):=NULL]
  new_valleys[,category:="UMR-DMV"]
  setkey(new_valleys,condition,start,end)
  
  # Add Valleys to UMRs and LMRs
  ov = foverlaps(seg_call,new_valleys,which=T,mult="all")
  seg_call=copy(seg_call[ov[is.na(yid),xid]])
  seg_call[,c("is.high","low.id","low.id.width","low.id.beta","segment_id"):=NULL]
  seg_call=rbind(seg_call,new_valleys,fill=T)
  seg_call[,estimate_id:=as.integer(condition)]
  setkey(seg_call,condition,start,end)
  seg_call[,segment_id:=as.numeric(seq_len(.N)),by=c("condition")]
  seg_call[,follows.large.gap:=start-shift(end)>pmd_valley_min_width,by=condition]
  seg_call[is.na(follows.large.gap),follows.large.gap:=F]
  
  
  if (split.pmds==T) {
    # split potential PMDs at UMRs and isolate from Valleys and large gaps
    umr_complement = seg_call[,.(estimate_id,start,end,category=ifelse((!is.na(category))&category%in%c("UMR","UMR-DMV"),"isolate","complement"),follows.large.gap)]
  } else {
    # Isolate PMDs from Valleys and large gaps
    umr_complement = seg_call[,.(estimate_id,start,end,category=ifelse((!is.na(category))&category=="UMR-DMV","isolate","complement"),follows.large.gap)]
  }
  umr_complement[,cat.chg:=shift(category)!=category,by=estimate_id]
  umr_complement[is.na(cat.chg),cat.chg:=F]
  umr_complement[,segment_id:=cumsum(cat.chg)+cumsum(follows.large.gap),by=estimate_id]
  umr_complement=umr_complement[,.(start=min(start),end=max(end)),by=c("estimate_id","segment_id","category")]
  umr_complement[,segment_id:=NULL]
  setkey(umr_complement,estimate_id,start,end)

  split_call = foverlaps(std_call,umr_complement)
  split_call = split_call[,.(start=pmax(start,i.start),end=pmin(end,i.end)),by=c("estimate_id","segment_id","i.start","i.end","fused_std","pseudoweight","category")]
  split_call[category!="complement",c("fused_std"):=list(0)]
  split_call[,c("i.start","i.end","segment_id","category"):=NULL]
  setkey(split_call,estimate_id,start,end)
  split_call[,segment_id:=seq_len(.N),by=estimate_id]

  
  # Gather individual methylation values for PMDs
  split_call <- split_call[start<=end]
  ov = foverlaps(mdata,split_call)
  ov[,c("start","end"):=list(min(pos),max(pos)+2),by=c("condition","estimate_id","segment_id")]
  pmd_call = ov[!is.na(segment_id),.(beta=mean(Nmeth/coverage),coverage=sum(coverage),pseudoweight=pseudoweight[1],
                                     num.cpgs=length(unique(pos)),betavals=list(Nmeth/coverage),coverages=list(coverage),fused_std=fused_std[1]),
                keyby=c("condition","estimate_id","segment_id","start","end")]
  pmd_call[,width:=end-start+1]

  # Call PMDs
  pmd_call[width>pmd_valley_min_width&width/(num.cpgs-1)<max_distance,is.admissible:=T]
  pmd_call[is.admissible==T&fused_std>=pmd_std_threshold&beta<pmd_max_beta,category:="PMD"]
  pmd_call[is.na(category),category:="other"]
  
  if ("segment_id"%in%names(pmd_call)) pmd_call[,segment_id:=NULL]
  if (merge.pmds==T) {
    # Merge PMDs
    pmd_call[,pmd.diff:=category!=shift(category),by=c("condition")]
    pmd_call[is.na(pmd.diff),pmd.diff:=TRUE]
    pmd_call[,follows.large.gap:=start-shift(end)>pmd_valley_min_width,by=condition]
    pmd_call[is.na(follows.large.gap),follows.large.gap:=F]
    pmd_call[,segment_id:=cumsum(pmd.diff)+cumsum(follows.large.gap),by=c("condition")]
  } else {
    pmd_call[,segment_id:=seq_len(.N),by=c("condition")]
  }
  
  # Annotate PMDs
  pmd_call = pmd_call[,.(start=min(start),end=max(end),betavals=list(unlist(betavals)),coverages=list(unlist(coverages)),
                         num.cpgs=sum(num.cpgs),fused_std=mean(fused_std)),
                      by=c("condition","category","segment_id")]
  pmd_call[,c("width","beta","std"):=list(end-start+1,mean(betavals[[1]]),sd(betavals[[1]])),
           by=c("condition","category","segment_id")]
  pmd_call[is.na(std),std:=0]
  pmd_call[,c("betavals","coverages"):=NULL]
  
  # Report PMD annotation
  setkeyv(seg_call,c("condition","start","end"))
  setkeyv(pmd_call,c("condition","start","end"))
  seg_call = foverlaps(seg_call,pmd_call[,.(condition,segment_id,start,end,category)])
  seg_call[,c("start","end","segment_id","estimate_id"):=NULL]
  setnames(seg_call,c("i.start","i.end","i.segment_id","category","i.category"),c("start","end","segment_id","PMD.status","category"))
  setkey(seg_call,condition,start,end)
  
  # Wrap up
  seg_call = unique(seg_call)
  seg_call[PMD.status=="PMD"&category=="LMR",category:=NA] #remove LMRs in PMDs
  seg_call[is.admissible==T,PMD.status:=ifelse(is.na(category)|category!="UMR",PMD.status,
                                               ifelse(shift(PMD.status)==shift(PMD.status,type="lead"),shift(PMD.status),"other")),by=condition] #set PMD status for UMRs based on neighboring segments
  if (drop==T) {
    seg_call[,c("is.V","is.flank","flank.dist","is.sep","fused_std"):=NULL]
    seg_call = seg_call[!is.na(category)]
    pmd_call = pmd_call[category=="PMD"]
  }
  # return
  segments=list(lmr_umr_valley=seg_call,pmd=pmd_call,ov=ov)
  return(segments)
}

#' Segment methylation data
#' 
#' @inheritParams segment_methylation_singlechr
#' @param ncores number of cores to parallelize on (default=1)
#'
#' @export
segment_methylation = function(data, ret, ncores=1,min_num_cpgs=4 , ...) {
  registerDoParallel(cores=ncores)
  segments = foreach (ch=data[,unique(chr)],.combine=function(x,y){list(lmr_umr_valley=rbind(x$lmr_umr_valley,y$lmr_umr_valley),pmd=rbind(x$pmd,y$pmd),ov=rbind(x$ov,y$ov))}) %dopar% {
    mdata=data[chr==ch]
    mret=list(estimates=ret$estimates[ret$estimates$chr==ch,], detection.type=ret$detection.type, ref.cond=ret$ref.cond)
    segments = MethyLasso:::segment_methylation_singlechr(mdata,mret, ...)
    segments$lmr_umr_valley[,chr:=ch]
    segments$pmd[,chr:=ch]
    segments$ov[,chr:=ch]
    segments
  }
  setcolorder(segments$lmr_umr_valley,c("condition","chr","start","end","width","segment_id","beta","coverage","pseudoweight","num.cpgs",
                                 "std","category","PMD.status"))
  setcolorder(segments$pmd,c("condition","chr","start","end","width","segment_id","beta","std","category"))
  return(segments)
}

#' call differences
#' 
#' single chromosome version
#' 
call_differences_singlechr = function(data, ret, min.diff=0.1,min.delta.diff=0.1, tol_val=0.01) {
  #call segments and bin original data along them
  if (ret$detection.type != "difference") stop("Please call call_differences on the result of difference_detection")
  mdata = copy(as.data.table(data))
  mdata[,coverage:=sum(coverage),by=c("condition","chr","pos")]
  mdata[,Nmeth:=sum(Nmeth),by=c("condition","chr","pos")]
  #mdata[,beta:=Nmeth/coverage]
  mdata[,beta:=mean(beta,na.rm=T),by=c("condition","chr","pos")]
  mdata=unique(subset(mdata,,c(condition,chr,pos,coverage,Nmeth,beta)))

  ref.cond=ret$ref
  diff.df=ret$estimates[ret$estimates$estimate_id > 1,] #drop common signal estimate
  if (is.character(ref.cond) || is.factor(ref.cond)) mdata[,is.ref:=condition==ref.cond] else mdata[,is.ref:=as.integer(condition)==ref.cond]
  mdata = foreach(cond=mdata[is.ref==F,unique(condition)],.combine=rbind) %do%
    rbind(mdata[condition==cond,.(is.ref,condition,pos,coverage,Nmeth,beta)],mdata[is.ref==T,.(is.ref,condition=cond,pos,coverage,Nmeth,beta)])
  mdata[,c("estimate_id","start","end"):=list(as.integer(factor(condition))+1,pos,pos+1)]
  setkey(mdata,estimate_id,start,end)
  
  diff_call = as.data.table(MethyLasso:::segment_methylation_helper(diff.df, tol_val = tol_val))
  stopifnot(diff_call[,uniqueN(estimate_id)] == mdata[,uniqueN(condition)])
  setkey(diff_call,estimate_id,start,end)

  diff_call <- diff_call[start<=end]
  ov = foverlaps(mdata,diff_call)
  max_distance=500

 # Split segments containing gaps
  ov[,is.oversize:=pos-shift(pos)>max_distance,keyby=c("condition","estimate_id","start","end")]
  ov[is.na(is.oversize),is.oversize:=F]
  ov[,segment_id:=segment_id+cumsum(is.oversize),keyby=c("condition","estimate_id")]

  # Shrink boundaries to span of CpGs
  ov[,c("start","end"):=list(min(pos),max(pos)+2),keyby=c("condition","estimate_id","segment_id")]
  ov[,width:=end-start+1]

  # Compute averages on each segment
  seg_call = ov[,.(beta=mean(i.beta,na.rm=T),coverage=sum(coverage),num.cpgs=length(unique(pos)),betavals=list(i.beta),coverages=list(coverage)),
                keyby=c("condition","estimate_id","segment_id","start","end","width","pseudoweight","is.ref")]
  seg_call[,is.ref:=ifelse(is.ref==T,".ref","")]
  seg_call = dcast(seg_call,condition+estimate_id+segment_id+start+end+width+pseudoweight~is.ref,
                   value.var = c("beta","coverage","num.cpgs","betavals","coverages"),sep = "")
  seg_call[is.na(coverage),c("coverage","num.cpgs"):=list(0,0)]
  seg_call[is.na(coverage.ref),c("coverage.ref","num.cpgs.ref"):=list(0,0)]
  seg_call[,diff:=beta-beta.ref]
  
  #group neighboring segments when their delta(diff) is small 
   seg_call[,group_change:=(!(abs(diff-shift(diff))<=min.delta.diff&abs(diff)>=min.diff&abs(shift(diff))>=min.diff&start-shift(end)<max_distance)),by=estimate_id]

  seg_call[is.na(group_change),group_change:=T]
  seg_call[,group_id:=cumsum(group_change),by=estimate_id]
  group_call=seg_call[,.(start=min(start),end=max(end),width=sum(width),coverage=sum(coverage),num.cpgs=sum(num.cpgs),
                         coverage.ref=sum(coverage.ref),pseudoweight=sum(pseudoweight),num.cpgs.ref=sum(num.cpgs.ref),
                         betavals=list(unlist(betavals)),betavals.ref=list(unlist(betavals.ref)),
			 coverages=list(unlist(coverages)),coverages.ref=list(unlist(coverages.ref))),
                      by=c("condition","estimate_id","group_id")]
  group_call[num.cpgs>0,c("beta","std"):=list(mean(betavals[[1]]),sd(betavals[[1]])),by=c("condition","estimate_id","group_id")]
  group_call[num.cpgs.ref>0,c("beta.ref","std.ref"):=list(mean(betavals.ref[[1]]),sd(betavals.ref[[1]])),by=c("condition","estimate_id","group_id")]
  group_call[,diff:=beta-beta.ref]
  group_call[,cohen.d:=diff/sqrt(((num.cpgs-1)*std**2+(num.cpgs.ref-1)*std.ref**2)/(num.cpgs+num.cpgs.ref-2))]
  setnames(group_call,"group_id","segment_id")
  #
  group_call[coverage>0 & coverage.ref>0, pval:=suppressWarnings(wilcox.test(betavals[[1]],betavals.ref[[1]],conf.int = F)$p.value),
           by=c("condition","estimate_id","segment_id","start","end","width")]
  group_call[is.na(pval),pval:=1] #ties with num.cpg=1 or coverage == 0 or coverage.ref == 0
  group_call[,c("betavals","betavals.ref","coverages","coverages.ref"):=NULL]

  cdata<- as.data.table(table(data$pos))
  cdata[,start:=as.numeric(V1)]
  cdata[,end:=as.numeric(start)+1]
  cdata[,score:=as.numeric(N)*100/length(levels(data$dataset))]
  setkey(cdata,start,end)
  ov3 <- foverlaps(group_call[,c("start","end")], cdata[,c("start","end","score")])
  setkeyv(ov3,c("i.start","i.end"))
  ov4 <- ov3[, .(coverage.score=mean(score)), by=c("i.start","i.end")]
  ov5 <- merge(group_call,ov4, by.x=c("start","end"),by.y=c("i.start","i.end"))
  group_call<- ov5
  
  #return only segments above threshold
  return(group_call)
}

#' call differences
#' 
#' @param data the input data frame
#' @param ret the output of difference_detection
#' @param min.diff minimum threshold for which a difference should be considered (default 0.1)
#' @param fdr.cutoff cutoff for FDR (default 0.05)
#' @param min.num.cpgs Require reported differences to have a certain minimum number of CpGs (default 4)
#' @param tol_val fused values below tol_val from each other are considered equal (default 0.01)
#' @param ncores number of cores to parallelize on (default=1)
#' 
#' @export
call_differences = function(data, ret, min.diff=0.1, pval.cutoff=0.05, fdr.cutoff=1, cov.score=70, min.num.cpgs=4, tol_val=0.01, ncores=1, verbose=TRUE) {
  registerDoParallel(cores=ncores)
  diff_call = foreach (ch=data[,unique(chr)],.combine=rbind) %dopar% {
  	if (verbose) cat("Processing chromosome:", ch, "\n")
    mdata=data[chr==ch]
    mret=list(estimates=ret$estimates[ret$estimates$chr==ch,], detection.type=ret$detection.type, ref=ret$ref)
    diff_call = MethyLasso:::call_differences_singlechr(mdata,mret,min.diff=min.diff,tol_val=tol_val)
    diff_call[,chr:=ch]
    diff_call[,ref:=ret$ref.cond]
    diff_call
    as.data.table(diff_call)
  }
  diff_call[,fdr:=p.adjust(pval,method='fdr')]
  setcolorder(diff_call,c("condition","ref","estimate_id","segment_id","chr","start","end","width","beta","std","coverage","num.cpgs",
                       "beta.ref","std.ref","coverage.ref","num.cpgs.ref","pseudoweight","diff","cohen.d","pval","fdr","coverage.score"))
  return(diff_call[(!is.na(diff))&abs(diff)>=min.diff&coverage.score>=cov.score&pval<=pval.cutoff&fdr<=fdr.cutoff&pmin(num.cpgs,num.cpgs.ref)>=min.num.cpgs])

}



#' Annotate regions with segmentation info
#'
#' @param regions the regions to annotate
#' @param segments the segments, as output by segment_methylation
#' @param min.overlap.pc differences which have no big overlap with any segment
#'   are classified as other (default 10 percent)
#'
annotate_regions = function(regions, segments, min.overlap.pc=10) {
  #compute PMD% for each difference
  pmd_call=segments$pmd[category=="PMD",.(condition,chr,start,end,category)]
  setkey(regions,condition,chr,start,end)
  setkey(pmd_call,condition,chr,start,end)
  ov2 = foverlaps(regions[,.(condition,chr,start,end,segment_id,diff.width=width-1)], pmd_call[,.(condition,chr,start,end,category)],nomatch=NULL,mult="all")
  ov2[,overlap:=pmin(end,i.end)-pmax(start,i.start)]
  ov2=ov2[,.(overlap=sum(overlap)),keyby=c("condition","chr","segment_id","diff.width")]
  ov2[,overlap.pc:=100*overlap/diff.width]
  ov2=ov2[,.(condition,chr,segment_id,pmd.pc=overlap.pc)]
  regions=ov2[regions,,on=c("condition","chr","segment_id")]
  regions[is.na(pmd.pc),pmd.pc:=0]
  #find overlaps with UMR/LMR/Valley
  sig_call=segments$lmr_umr_valley[!is.na(category),.(condition,chr,start,end,category)]
  setkey(sig_call,condition,chr,start,end)
  ov = foverlaps(regions[,.(condition,chr,start,end,segment_id,diff.width=width-1,pmd.pc)],
                 sig_call[,.(condition,chr,start,end,category)],
                 nomatch=NA,mult="all")
  #count number of UMRs/LMRs which overlap
  ov[,c("num.UMR","num.LMR"):=list(sum(grepl("UMR",category)),sum(grepl("LMR",category))),by=c("condition","chr","segment_id")]
  #count which type of segment overlaps most with that difference, and keep largest overlap
  #annotate non-overlapping segments as "other" or "PMD" based on pmd.pc
  ov[,overlap:=pmin(end,i.end)-pmax(start,i.start)]
  ov[is.na(overlap),overlap:=0]
  ov[is.na(category),category:=ifelse(pmd.pc<50,"other","PMD")]
  ov=ov[,.(overlap=sum(overlap)),keyby=c("condition","chr","segment_id","diff.width","num.UMR","num.LMR","category","pmd.pc")]
  ov[,overlap.pc:=100*overlap/diff.width]
  setkeyv(ov,c("condition","chr","segment_id","overlap.pc"))
  ov[,overlap.id:=seq.int(.N,1),by=c("condition","chr","segment_id")]
  ov=ov[overlap.id==1]
  #differences which have no big overlap with any segment are classified as "other" or "PMD"
  ov[overlap.pc<min.overlap.pc,category:=ifelse(pmd.pc<50,"other","PMD")]
  ov=ov[,.(condition,chr,segment_id,category,num.UMR,num.LMR)]
  return(ov)
}

#' Annotate differences with segmentation info
#'
#' @param diff_call the differences data frame, as returned by call_differences
#' @param segments the segments, as output by segment_methylation
#' @param min.overlap.pc differences which have no big overlap with any segment
#'   are classified as "other" (default 10 percent)
#'   
#' @export
#'
annotate_differences = function(diff_call, segments, min.overlap.pc=10) {
  # Annotate with target segments
  ov.target = MethyLasso:::annotate_regions(diff_call[,.(condition,chr,start,end,segment_id,width)],segments,min.overlap.pc=min.overlap.pc)
  ov.target = diff_call[,.SD[1,.(ref)],by=condition][ov.target,,on="condition"]

  # Annotate with reference segments
  ov.ref = MethyLasso:::annotate_regions(diff_call[,.(condition=ref,chr,start,end,segment_id,width)],segments,min.overlap.pc=min.overlap.pc)
  setnames(ov.ref,"condition","ref")
  ov.ref = diff_call[,.SD[1,.(condition)],by=ref][ov.ref,,on="ref"]
  
  # Merge
  ov = merge(ov.target,ov.ref,all=T,by=c("condition","ref","chr","segment_id"),suffixes = c("",".ref"))
  ov[,switch:=ifelse(category.ref==category, paste("within",category.ref),
                     ifelse(category.ref=="other",paste("to",category),
                            ifelse(category=="other",paste("from",category.ref),
                                   paste(category.ref,"->",category))))]
  ov = ov[,.(condition,ref,chr,segment_id,category,category.ref,switch,delta.UMR=num.UMR-num.UMR.ref,delta.LMR=num.LMR-num.LMR.ref)]
  differences=ov[diff_call,,on=c("condition","ref","chr","segment_id")]
  setcolorder(differences,c("condition","ref","estimate_id","segment_id","chr","start","end","width","beta","std","coverage","num.cpgs",
                           "beta.ref","std.ref","coverage.ref","num.cpgs.ref","pseudoweight","diff","cohen.d","pval","fdr","category","category.ref",
                          "switch","delta.UMR","delta.LMR"))
  # Return
  return(differences)
}

