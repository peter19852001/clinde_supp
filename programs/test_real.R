###
## To test on real data from TimeDelayARACNE
## and compare the predicted structure and the true structure
## Since the delay for the links are not known, and for some, the effect (+/-)
## is not known, we compare only the performance on the links and effect.

# assume test.n, test.p, test.m, test.max.n are defined

source("sim.R");
source("infer.R");
#
cal_f_measure <- function(r,p) {2*p*r/(p+r)}
#
linear.interpolate.exp <- function(n.p,r) {
  # r is n by g matrix, where n is number of time points, g is number of genes
  # n.p is the ratio of oversampling
  # returns a matrix of n.p*n by r, where the columns are separately oversampled linearly
  n <- nrow(r);
  z <- matrix(rep(0,n.p*n*ncol(r)), nrow=n.p*n, ncol=ncol(r));
  for(i in 1:ncol(r)) {
    tmp <- approx(x=1:n, y=r[,i], n=n*n.p);
    z[,i] <- tmp$y;
  }
  z
}
spline.interpolate.exp <- function(n.p,r) {
  # r is n by g matrix, where n is number of time points, g is number of genes
  # n.p is the ratio of oversampling
  # returns a matrix of n.p*n by r, where the columns are separately oversampled using natural spline
  n <- nrow(r);
  z <- matrix(rep(0,n.p*n*ncol(r)), nrow=n.p*n, ncol=ncol(r));
  for(i in 1:ncol(r)) {
    tmp <- spline(x=1:n, y=r[,i], n=n*n.p, method="natural");
    z[,i] <- tmp$y;
  }
  colnames(z) <- colnames(r);
  z
}

interpolate.exp.sec <- function(n.p, r, f, segs=NULL) {
  # r is n by g matrix, where n is number of time points, g is number of genes
  # n.p is the ratio of oversampling
  # f is the function to sample each section, e.g. spline.interpolate.exp
  # segs records the lengths of each section
  if(is.null(segs)) {
    f(n.p, r)
  } else {
    z <- f(n.p,  r[1:segs[1],]);
    if(length(segs) > 1) {
      L <- segs[1];
      for(i in 2:length(segs)) {
        z <- rbind(z, f(n.p, r[(L+1):(L+segs[i]),]));
        L <- L + segs[i];
      }
    }
    z
  }
}
#
test1 <- function(data_name, true_grn,real_data, p1,p2,max_n,max.delay, take.log, m.interpolate, n.interpolate, method=c("pcor","mi"), pruning=c("all","common"), segs=NULL, one.delay=FALSE, no.dup=FALSE) {
  # true_grn is a matrix of n by n, where n is the number of genes
  #    true_grn[i,j] entry is the causal effect of gene i to gene j
  # use the expression data in real_data to infer a gene network
  # the expression data is in suitable format for infer_grn1.
  # if take.log is true, then take log for the expression before interpolation and further processing
  # m.interpolate is either "none", "linear" or "spline", determining the method of interpolation
  # n.interpolate is the how many times to oversample by either linear or spline interpolation (natural spline)
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  if(method == "pcor") {u.test.f <- test.corr; c.test.f <- test.p.corr;}
  if(method == "mi") {u.test.f <- mutual_info; c.test.f <- c_mutual_info;}
  #
  tmp.sr <- if(take.log) log(0.00001 + real_data) else real_data; # try taking log of the real expression data
  if(m.interpolate == "linear") {
    sr <- interpolate.exp.sec(n.interpolate, tmp.sr, linear.interpolate.exp, segs);
  } else if(m.interpolate == "spline") {
    sr <- interpolate.exp.sec(n.interpolate, tmp.sr, spline.interpolate.exp, segs);
  } else {
    sr <- tmp.sr;
  }
  n.timepoints <- nrow(tmp.sr);
  n.interpolate.points <- n.timepoints * n.interpolate;
  # printing some information
  cat("\n## data=",data_name,"\n");
  cat("\n## p1=",p1,"\n");
  cat("\n## p2=",p2,"\n");
  cat("\n## max.n=",max_n,"\n");
  cat("\n## max.delay=",max.delay,"\n");
  cat("\n## max.delay * n.interpolate=",max.delay*n.interpolate,"\n");
  cat("\n## n.points=",n.timepoints,"\n");
  cat("\n## n.points (after interpolation)=",n.interpolate.points,"\n");
  cat("\n## take.log=",take.log,"\n");
  cat("\n## m.interpolate=",m.interpolate,"\n");
  cat("\n## n.interpolate=",n.interpolate,"\n");
  cat("\n## method=",method,"\n");
  cat("\n## pruning=",pruning,"\n");
  #
  pn1 <- infer_grn1(sr,p1,n.interpolate*max.delay, u.test.f, n.interpolate*segs); # max.delay is counted in time steps, adjusted for interpolation
  if(one.delay) {pn1 <- remove.dup.links(pn1);}
  pn2 <- filter_pcor(sr,p2,pn1,max_n, c.test.f, pruning, n.interpolate*segs);
  if(no.dup) {pn2 <- remove.dup.links(pn2);}
  # stage 1
  cat("== Stage 1\n");
  print(pn1);
  r1 <- compare_grn(true_grn, true_grn, pn1); # for the delay, just put dummy there, ignore it
  cat("== Stage 1 Performance\n");
  print(r1);
  # stage 2
  cat("== Stage 2\n");
  print(pn2);
  r2 <- compare_grn(true_grn, true_grn, pn2); # for the delay, just put dummy there, ignore it
  cat("== Stage 2 Performance\n");
  print(r2);
  # print one line comma separated summary for the various metrics, in particular ignore the metrices on delay
  cat("\n####,",
      data_name,",",p1,",",p2,",",max_n,",",max.delay,",",n.timepoints,",",
      n.interpolate.points,",",take.log,",",m.interpolate,",",n.interpolate,",",method,",",pruning,",",
      r1$links.recall,",",r1$links.precision,",",cal_f_measure(r1$links.recall,r1$links.precision),",",
      r2$links.recall,",",r2$links.precision,",",cal_f_measure(r2$links.recall,r2$links.precision),",",
      r1$effects.recall,",",r1$effects.precision,",",cal_f_measure(r1$effects.recall,r1$effects.precision),",",
      r2$effects.recall,",",r2$effects.precision,",",cal_f_measure(r2$effects.recall,r2$effects.precision),"\n",
      sep="");
  #
  list(stage1.grn=pn1, stage2.grn=pn2, stage1.res=r1, stage2.res=r2)
}

grn.to.matrix <- function(g,n) {
  # g is the network as returned by infer_grn1 or infer_grn.
  # n is the number of genes
  # represent g in the form of a matrix, where the (i,j)th entry is test.value of gi to gj.
  # if there are multiple entries for an pair (i,j), take the first one.
  r <- matrix(rep(0,n*n), nrow=n, ncol=n);
  for(i in 1:nrow(g)) {
    x <- g$from[i];
    y <- g$to[i];
    if(r[x,y] == 0) {
      r[x,y] <- g$test.value[i];
    }
  }
  r
}

standardize.exp <- function(r) {
  # each column of r is expression, standardize each column, so that each one has variance 1 by dividing by the s.d.
  for(i in 1:ncol(r)) {
    x <- sd(r[,i]);
    if(x > 0) {r[,i] <- r[,i]/x;}
  }
  r
}
# load the data, which contains the different datasets we need
load("real_data.RData");

load("exp_data/kegg_meiosis_exp.RData");
load("kegg_meiosis_sn.RData");
# a call should be made here.
# e.g. test1("yeast",yeast_true,yeast_exp, 0.0001,0.001,4,3)
# test1("yeast",yeast_true,yeast_exp, 0.9, 0.99, 2,4, FALSE,"linear",1,"pcor")
# test1("yeast",yeast_true,yeast_exp, 0.3, 0.3, 2,4, FALSE,"linear",1,"mi")
# test1("SOS",sos_true,sos_exp, 0.9, 0.5, 2,2, FALSE,"linear",1,"pcor") 
# test1("SOS",sos_true,sos_exp, 0.3, 0.3, 2,4, FALSE,"linear",1,"mi") 
# test1("IRMA",irma_true,irma_exp, 0.9, 0.99, 2,4, FALSE,"linear",1,"pcor") 
# test1("SOS",sos_true,sos_exp, 0.3, 0.3, 2,4, FALSE,"linear",1,"mi")

# test1("kegg_meiosis", kegg.meiosis.true, kegg.meiosis.occ.exp, 3,3,4,4, FALSE, "linear",1, "pcor", "all", kegg.meiosis.i.seg, FALSE,TRUE);
#test1("kegg_meiosis_sn1", kegg.meiosis.true.sn1, kegg.meiosis.exp.i.sn1, 1.5,1.5,4,4, FALSE, "linear",1, "pcor", "all", kegg.meiosis.i.seg, FALSE, TRUE);
#test1("kegg_meiosis_sn1", kegg.meiosis.true.sn1, kegg.meiosis.exp.i.sn1.small, 1.5,1.5,4,4, FALSE, "linear",1, "pcor", "all", c(40), FALSE, TRUE);
