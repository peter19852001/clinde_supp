#
# to randomly generate (synthetic) gene network in the form of matrices,
# where each link a_ij has a number representing the
# effect of gene i on gene j: +ve for activation, -ve for inhibition.
# Associated with each a_ij =/= 0 is t_ij > 0, which represents the time delay
# of the effect of gene i on gene j.

gen_grn <- function(n,is_acyclic,is_self, max.parents) {
  # n is the number of gene
  # is_acyclic is true iff only acyclic network is to be generated
  # is_self is true iff when is_acyclic is false, to include self loops
  # p_link is the probability of a link
  # Returns a list of two matrices of n by n, one is the links, the other
  # is the delays.
  L <- n*n;
  r <- rep(0,L);
  # limits the number of non-zero entries in each column
  for(j in 1:n) {
    ps <- sample(1:n,max.parents);
    for(i in 1:length(ps)) {
      si <- (ps[i]-1)*n + j;
      r[si] <- 1;
    }
  }
  #
  if(is_acyclic) {
    # make it upper triangular
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        r[(i-1)*n + j] <- 0;
      }
    }
  }
  if(!is_self) {
    # make the diagonal zero
    for(i in 1:n) {r[(i-1)*n + i] <- 0;}
  }
  delay <- r*runif(L,0,1); # simulate less delays
  r <- r*(1-2*rbinom(L,1,0.5))*runif(L,0.5,1.5);
  # returns r and delay
  list(links=matrix(data=r, nrow=n, ncol=n, byrow=TRUE),
       delays=matrix(data=delay, nrow=n, ncol=n, byrow=TRUE))
}

permute_grn <- function(links, delays) {
  # to randomly permute the genes so that the position do not give any advantage
  x <- sample(1:nrow(links));
  list(links=links[x,x], delays=delays[x,x])
}

sim_grn <- function(links, delays, dt, N) {
  # links is a matrix encoding the links of a gene network (direction and magnitude)
  # delays contains the time delay for each link
  # dt is the time step to simulate, and to simulate N steps
  # Returns a matrix a matrix of N by n, where n is the number of genes in the gene network.
  # The network is assumed to start with zero expression for each gene
  T <- ceiling(delays/dt); # turn the delays into steps
  ng <- nrow(links);
  r <- matrix(data=rep(0,ng*N),nrow=N,ncol=ng);
  for(i in 1:N) {
    for(j in 1:ng) {
      x <- rnorm(1,0,0.01);
      for(k in 1:ng) {
        if(T[k,j] != 0) {
          x <- x + (if(T[k,j] < i) r[i-T[k,j],k] else 0)*links[k,j];
        }
      }
      r[i,j] <- x;
    }
  }
  r
}

plot_exp <- function(r) {
  # r is a n by g matrix, where n is the number of time points,
  # g is the number of genes
  # The values are the expression of the genes at the different time points
  n <- nrow(r);
  g <- ncol(r);
  legend.r <- if(is.null(colnames(r))) rep("",g) else colnames(r);
  for(i in 1:g) {legend.r[i] <- paste("i=",i,legend.r[i]);}
  plot(x=c(1,n),y=range(r), type="n", main="Expressions",xlab="Time",ylab="Expression");
  for(i in 1:g) {
    lines(x=1:n, y=r[,i], type="b",col=i,pch=i);
  }
  legend(x="topright",legend=legend.r,pch=1:g,col=1:g);
}


compare_matrices <- function(a,b) {
  # a and b are two matrices of the same size
  # report the entries for which a and b are different
  for(i in 1:nrow(a)) {
    for(j in 1:ncol(b)) {
      if(a[i,j] != b[i,j]) {cat("i: ",i, "\tj: ",j,"\ta: ",a[i,j], "\tb: ",b[i,j],"\n");}
    }
  }
}

compare_grn_old <- function(tlinks,tdelays, elinks,edelays) {
  # tlinks and tdelays are the true links and delays, respectively
  # elinks and edelays are the estimated links and delays
  # all four are square matrices of the size n by n, where n is the number of genes
  # Returns a list of:
  #  delays.right: among the positives (either true or predicted), the number of delays correctly estimated
  #  delays.wrong: among the positives (either true or predicted), the number of delays incorrectly estimated
  #  delays.sse: sum of squared errors of the delays in the TP
  #  links: the contingency table of the presence and absence of links (considering also the direction)
  #  effects: the contingency table of the signs (1 for +, 0 for no links, -1 for -) of the links
  #    with direction considered,
  #    where the rows are true signs, and columns are predicted signs,
  # delays
  TPs <- (sign(tlinks) != 0) | (sign(elinks) != 0);
  dr <- sum(TPs & (tdelays == edelays));
  dw <- sum(TPs & (tdelays != edelays));
  dsse <- (tdelays - edelays)[TPs];
  dsse <- sum(dsse*dsse);
  # links
  links.tab <- table(sign(tlinks)!=0,sign(elinks)!=0, dnn=c("True Links","Predicted Links"));
  # effects
  effects.tab <- table(sign(tlinks),sign(elinks), dnn=c("True Effect","Predicted Effect"));
  #
  list(delays.right=dr,delays.wrong=dw, delays.sse=dsse, links=links.tab, effects=effects.tab)
}

compare_grn <- function(tlinks,tdelays, grn) {
  # tlinks and tdelays are the true links and delays, respectively,
  #  both are square matrices of the size n by n, where n is the number of genes
  # grn is the result as returned by infer_grn(), and contains the predicted edges
  # Returns a list of:
  #  links.recall: the recall of the links, consider the direction (x -> y or y -> x) of edge,
  #    but not the effect (+/-), and disregard the delay
  #  links.precision: similar to links.recall, but for the precision of the links
  #  effects.recall: consider the direction and sign of the effect. The recall of the effects
  #  effects.precision: similar to effects.recall, but for precision
  #  delays.recall: among the delays for true links, how many are correctly predicted (same value)
  #  delays.precision: among the predicted delays, how many are correct (same value)
  ng <- nrow(tlinks);

  n.links <- sum(sign(tlinks)!=0); # same as number of effects, and number of delays
  n.nlinks <- sum(sign(tlinks)==0); # number of non-links
  n.p.links <- nrow(grn); # number of predicted links, same as number of predicted effects

  n.r.links <- 0; # number of true links correctly predicted
  n.c.links <- 0; # number of correct prediction in the links, multiple prediction for the same link (possibly with different delays) are all counted
  n.n.links <- 0;
  #
  n.r.effects <- 0; # number of true effects correctly predicted
  n.c.effects <- 0; # number of correct prediction of the effects, with proper multiple counts
  n.r.delays <- 0; # number of true delays correctly predicted
  n.c.delays <- 0; # number of correct prediction of the delays, with proper multiple counts
  # go through the true links
  if(nrow(grn) > 0) {
    for(i in 1:ng) {
      for(j in 1:ng) {
        s <- sign(tlinks[i,j]);
        if(s != 0) {
          d <- (grn$from==i) & (grn$to==j);
          if(sum(d) > 0) {n.r.links <- n.r.links + 1;}
          if(sum(d & (grn$delay == tdelays[i,j])) > 0) {n.r.delays <- n.r.delays + 1;}
          if(s < 0) {
            if(sum(d & (grn$test.value < 0)) > 0) {n.r.effects <- n.r.effects + 1;}
          } else {
            if(sum(d & (grn$test.value > 0)) > 0) {n.r.effects <- n.r.effects + 1;}
          }
        } else {
          d <- (grn$from==i) & (grn$to==j);
          if(sum(d) == 0) {n.n.links <- n.n.links + 1;}
        }
      }
    }
  }
  # go through the predictions
  if(nrow(grn) > 0) {
    for(i in 1:nrow(grn)) {
      x <- grn$from[i];
      y <- grn$to[i];
      if(tlinks[x,y] != 0) {n.c.links <- n.c.links + 1;}
      if(grn$test.value[i] < 0) {
        if(tlinks[x,y] < 0) {n.c.effects <- n.c.effects + 1;}
      } else {
        if(tlinks[x,y] > 0) {n.c.effects <- n.c.effects + 1;}
      }
      if(tdelays[x,y] == grn$delay[i]) {n.c.delays <- n.c.delays + 1;}
    }
  }
  # done
  list(links.recall=n.r.links/n.links, links.precision=(if(n.p.links<=0) 0 else (n.c.links/n.p.links)),
       links.specificity=n.n.links/n.nlinks,
       effects.recall=n.r.effects/n.links, effects.precision=(if(n.p.links<=0) 0 else (n.c.effects/n.p.links)),
       delays.recall=n.r.delays/n.links, delays.precision=(if(n.p.links<=0) 0 else (n.c.delays/n.p.links)))
}

### test
#tmp <- gen_grn(20,TRUE,FALSE,0.2);
#tmp$delays <- ceiling(tmp$delays*10);
#tmpr <- sim_grn(tmp$links,tmp$delays,0.1,1000);
#plot_exp(tmpr);
#tmpz1 <- infer_grn1(tmpr,0.0001,100);
#tmpz2 <- infer_grn(tmpr,0.001,100);
#tmpz3 <- infer_grn(tmpr,0.01,100);
#tmpz <- infer_grn(tmpr,0.0001,100);

#compare_matrices(ceiling(tmp$delays*10), tmpz1$delays);
#compare_matrices(ceiling(tmp$delays*10), tmpz2$delays);
#compare_matrices(ceiling(tmp$delays*10), tmpz3$delays);
#compare_matrices(ceiling(tmp$delays*10), tmpz$delays);

#r1 <- compare_grn(tmp$links,tmp$delays, tmpz1$links,tmpz1$delays);
#r1
#r2 <- compare_grn(tmp$links,tmp$delays, tmpz2$links,tmpz2$delays);
#r2
#r3 <- compare_grn(tmp$links,tmp$delays, tmpz3$links,tmpz3$delays);
#r3
#r4 <- compare_grn(tmp$links,tmp$delays, tmpz$links,tmpz$delays);
#r4

#tmp
#tmpz

#ceiling(tmp$delays*10)
#tmpz1$delays

#ceiling(tmp$delays*10)
#tmpz2$delays

#ceiling(tmp$delays*10)
#tmpz3$delays

#ceiling(tmp$delays*10)
#tmpz$delays

##
