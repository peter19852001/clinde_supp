###
## To test the sim.R, i.e. randomly generating gene network,
## then generate the expression, then infer the structure,
## and compare the predicted structure and the true structure

# assume test.n, test.p, test.m, test.max.n are defined

source("sim.R");
source("infer.R");

# now for acyclic, no self loop case
test1 <- function(n,st1,st2,max_n, max.parents, max.delay, n.points, is.acyclic, is.self, method=c("pcor","mi"), pruning=c("all","common")) {
  # simulate and test for one random gene network of n genes, with p as the p.value threshold
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  if(method == "pcor") {u.test.f <- test.corr; c.test.f <- test.p.corr;}
  if(method == "mi") {u.test.f <- mutual_info; c.test.f <- c_mutual_info;}
  #
  tn <- gen_grn(n,is.acyclic,is.self, max.parents);
  tn <- permute_grn(tn$links, tn$delays);
  sr <- sim_grn(tn$links, tn$delays, 0.1, n.points);
  tn$delays <- ceiling(tn$delays*10);

  pn1 <- infer_grn1(sr,st1,max.delay, u.test.f); # max.delay is counted in time steps
  pn2 <- filter_pcor(sr,st2,pn1,max_n, c.test.f, pruning);
  # stage 1
  r1 <- compare_grn(tn$links, tn$delays, pn1);
  cat("== Stage 1\n");
  print(r1);
  # stage 2
  r2 <- compare_grn(tn$links, tn$delays, pn2);
  cat("== Stage 2\n");
  print(r2);
  list(stage1=r1, stage2=r2)
}

test.many <- function(n,st1,st2,m, max_n, max.parents,max.delay, n.points,is.acyclic,is.self, method=c("pcor","mi"), pruning=c("all","common")) {
  # call test1 m times
  # store the performance of each trial, and finally print the summary of the performances of all the trials
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  # for stage 1
  d1.recalls <- rep(0,m);
  d1.precisions <- rep(0,m);
  d1.f.measures <- rep(0,m);
  e1.recalls <- rep(0,m);
  e1.precisions <- rep(0,m);
  e1.f.measures <- rep(0,m);
  l1.recalls <- rep(0,m);
  l1.precisions <- rep(0,m);
  l1.f.measures <- rep(0,m);
  # for stage 2
  d2.recalls <- rep(0,m);
  d2.precisions <- rep(0,m);
  d2.f.measures <- rep(0,m);
  e2.recalls <- rep(0,m);
  e2.precisions <- rep(0,m);
  e2.f.measures <- rep(0,m);
  l2.recalls <- rep(0,m);
  l2.precisions <- rep(0,m);
  l2.f.measures <- rep(0,m);
  #
  for(i in 1:m) {
    cat("\n# Trial ",i,"\n");
    x <- test1(n,st1,st2,max_n, max.parents,max.delay, n.points,is.acyclic, is.self, method, pruning);
    s1 <- x$stage1;
    s2 <- x$stage2;
    # for stage 1
    d1.recalls[i] <- s1$delays.recall;
    d1.precisions[i] <- s1$delays.precision;
    d1.f.measures[i] <- 2*(s1$delays.recall)*(s1$delays.precision)/(s1$delays.recall + s1$delays.precision);
    e1.recalls[i] <- s1$effects.recall;
    e1.precisions[i] <- s1$effects.precision;
    e1.f.measures[i] <- 2*(s1$effects.recall)*(s1$effects.precision)/(s1$effects.recall + s1$effects.precision);
    l1.recalls[i] <- s1$links.recall;
    l1.precisions[i] <- s1$links.precision;
    l1.f.measures[i] <- 2*(s1$links.recall)*(s1$links.precision)/(s1$links.recall + s1$links.precision);
    # for stage 2
    d2.recalls[i] <- s2$delays.recall;
    d2.precisions[i] <- s2$delays.precision;
    d2.f.measures[i] <- 2*(s2$delays.recall)*(s2$delays.precision)/(s2$delays.recall + s2$delays.precision);
    e2.recalls[i] <- s2$effects.recall;
    e2.precisions[i] <- s2$effects.precision;
    e2.f.measures[i] <- 2*(s2$effects.recall)*(s2$effects.precision)/(s2$effects.recall + s2$effects.precision);
    l2.recalls[i] <- s2$links.recall;
    l2.precisions[i] <- s2$links.precision;
    l2.f.measures[i] <- 2*(s2$links.recall)*(s2$links.precision)/(s2$links.recall + s2$links.precision);
  }
  # For stage 1
  cat("\n====================\nStage 1\n====================\n");
  #
  cat("\nDelays Recall\n");
  print(summary(d1.recalls));
  cat("\nDelays Precision\n");
  print(summary(d1.precisions));
  cat("\nDelays F-measure\n");
  print(summary(d1.f.measures));
  #
  cat("\nLinks Recall\n");
  print(summary(l1.recalls));
  cat("\nLinks Precision\n");
  print(summary(l1.precisions));
  cat("\nLinks F-measure\n");
  print(summary(l1.f.measures));
  #
  cat("\nEffects Recall\n");
  print(summary(e1.recalls));
  cat("\nEffects Precision\n");
  print(summary(e1.precisions));
  cat("\nEffects F-measure\n");
  print(summary(e1.f.measures));
  # Stage 2
  cat("\n====================\nStage 2\n====================\n");
  cat("\nDelays Recall\n");
  print(summary(d2.recalls));
  cat("\nDelays Precision\n");
  print(summary(d2.precisions));
  cat("\nDelays F-measure\n");
  print(summary(d2.f.measures));
  #
  cat("\nLinks Recall\n");
  print(summary(l2.recalls));
  cat("\nLinks Precision\n");
  print(summary(l2.precisions));
  cat("\nLinks F-measure\n");
  print(summary(l2.f.measures));
  #
  cat("\nEffects Recall\n");
  print(summary(e2.recalls));
  cat("\nEffects Precision\n");
  print(summary(e2.precisions));
  cat("\nEffects F-measure\n");
  print(summary(e2.f.measures));
  # print one line comma separated summary for the various medians
  cat("\n####,",
      n,",",st1,",",st2,",",m,",", max_n,",",max.parents,",",max.delay,",",
      n.points,",", is.acyclic,",",is.self,",",method,",",pruning,",",
      median(e1.recalls, na.rm=TRUE),",",median(e1.precisions, na.rm=TRUE),",",median(e1.f.measures, na.rm=TRUE),",",
      median(e2.recalls, na.rm=TRUE),",",median(e2.precisions, na.rm=TRUE),",",median(e2.f.measures, na.rm=TRUE),",",
      median(l1.recalls, na.rm=TRUE),",",median(l1.precisions, na.rm=TRUE),",",median(l1.f.measures, na.rm=TRUE),",",
      median(l2.recalls, na.rm=TRUE),",",median(l2.precisions, na.rm=TRUE),",",median(l2.f.measures, na.rm=TRUE),",",
      median(d1.recalls, na.rm=TRUE),",",median(d1.precisions, na.rm=TRUE),",",median(d1.f.measures, na.rm=TRUE),",",
      median(d2.recalls, na.rm=TRUE),",",median(d2.precisions, na.rm=TRUE),",",median(d2.f.measures, na.rm=TRUE),"\n");
}

# assume test.n, test.p1, test.p2, test.m, test.max.n, test.max.parents, test.max.delay, test.n.points,
#   test.is.acyclic, test.is.self, test.method are defined
cat("\n## n=",test.n,"\n");
cat("\n## st1=",test.p1,"\n");
cat("\n## st2=",test.p2,"\n");
cat("\n## m=",test.m,"\n");
cat("\n## max.n=",test.max.n,"\n");
cat("\n## max.parents=",test.max.parents,"\n");
cat("\n## max.delay=",test.max.delay,"\n");
cat("\n## n.points=",test.n.points,"\n");
cat("\n## is.acyclic=",test.is.acyclic,"\n");
cat("\n## is.self=",test.is.self,"\n");
cat("\n## method=",test.method,"\n");
cat("\n## pruning=",test.pruning,"\n");
test.many(test.n,test.p1,test.p2,test.m, test.max.n, test.max.parents, test.max.delay, test.n.points, test.is.acyclic, test.is.self,test.method, test.pruning);
