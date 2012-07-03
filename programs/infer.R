
init.change.exp <- function(x, threshold.up, threshold.down) {
  # x is a vector of values, determine the minimum position i such that
  # x[1]/x[i] >= threshold.up or x[1]/x[i] <= threshold.down
  # Returns the minimum position i if one exists, otherwise, returns the last one
  for(i in 2:length(x)) {
    tmp <- x[1]/x[i];
    if(!is.na(tmp) && ((tmp >= threshold.up) || (tmp <= threshold.down))) {return(i)}
  }
  length(x)
}

test.link <- function(x,y,lags,st,test.f,segs=NULL) {
  # x and y are two vector (should have the same length), representing two time-series
  # lags is the maximum allowed time lags
  # st is the score threshold, below which the link is not significant
  # test.f is a function accepting two vectors (of the same length), and
  #  returns a list of (score,test.value), where the score is higher the better,
  #  test.value is whatever value representing the strenth of the test (maybe test statistic, or p-value).
  # segs is a vector of the length of the different segments in x and y, test.f is applied to each
  # the segments are shifted separately, then pieced together to calculate the score as a whole.
  # calculate the most significant score of the shifted x (x precede y).
  # Returns a list of lists of (score,test.value,delay), where test.value is the significant value given by test.f
  #  and delay is a non-negative integer: y[delay+i] <-> x[i].
  #  The returned list contains all the significant lags. An empty list is returned if there is
  #  no significant lags.
  if(is.null(segs)) {segs <- c(min(length(x),length(y)))}
  rL <- list();
  ri <- 1;
  for(i in 1:min(lags,min(segs)-5)) { # the delays should not go across segments
    #
    tmp.L <- 0;
    tmp.len <- min(length(x),length(y)) - (i-1)*length(segs);
    tmp.x <- rep(0,tmp.len);
    tmp.y <- rep(0,tmp.len);
    tmp.p <- 1;
    # now collect the different shifted segments
    for(j in 1:length(segs)) {
      L <- segs[j] - i + 1;
      tmp.x[tmp.p:(tmp.p+L-1)] <- x[(tmp.L+1):(tmp.L+L)];
      tmp.y[tmp.p:(tmp.p+L-1)] <- y[(tmp.L+i):(tmp.L+i+L-1)];
      #
      tmp.L <- tmp.L + segs[j];
      tmp.p <- tmp.p + L;
    }
    # can caluclate the score
    r <- test.f(tmp.y, tmp.x);
    if(is.numeric(r$score) && (r$score >= st)) {
      rL[[ri]] <- list(score=r$score, test.value=r$test.value, delay=i-1);
      ri <- ri+1;
    }
  }
  rL
}

infer_grn1 <- function(r,st,lags, test.f, segs=NULL) {
  # Infers the time lag, and pairwise correlation.
  # r is a n by g matrix, where n is the number of time points,
  # g is the number of genes
  # The values are the expression of the genes at the different time points.
  # st is the threshold for score, above which the link is significant
  # test.f is a function accepting two vectors (of the same length), and
  #  returns a list of lists of (score,test.value,delay), where the score is higher the better,
  #  test.value is whatever value representing the strenth of the test (maybe test statistic, or p-value).
  # lags is the maximum allowed time lag in inferring the correlation
  # Use these expressions to calculate the n by n (indirect links) between the
  #  genes and the associated time delay, when the score of the particular
  #  delayed link is >= st.
  # Returns a data frame with columns: from, to, score, test.value, delay,
  #  where the rows are the predicted edges. The edges may have different entries for the
  #  same pair of vertices, but with different delay values

  n <- nrow(r);
  g <- ncol(r);
  tmp <- list(); # used to temporarily store the significant delays
  tmpi <- 1;
  n.entries <- 0; # the number of entries in the final data frame
  for(i in 1:g) {
    for(j in 1:g) {
      if(i != j) { # currently no self loops considered
        z <- test.link(r[,i],r[,j],lags,st, test.f, segs);
        if(length(z) > 0) { # has entries, add to tmp
          tmp[[tmpi]] <- list(fi=i,ti=j,delays=z);
          tmpi <- tmpi + 1;
          n.entries <- n.entries + length(z);
        }
      }
    }
  }
  # now turn the entries into a data frame
  v.from <- rep(0,n.entries);
  v.to <- rep(0,n.entries);
  v.score <- rep(0,n.entries);
  v.test.value <- rep(0,n.entries);
  v.delay <- rep(0,n.entries);
  # get the values in tmp into the above vectors
  vi <- 1;
  if(length(tmp)>0) {
    for(i in 1:length(tmp)) {
      z <- tmp[[i]];
      for(j in 1:length(z$delays)) {
        v.from[vi] <- z$fi;
        v.to[vi] <- z$ti;
        x <- z$delays[[j]];
        v.score[vi] <- x$score;
        v.test.value[vi] <- x$test.value;
        v.delay[vi] <- x$delay;
        vi <- vi + 1;
      }
    }
  }
  # done
  s <- sort(v.score,index.return=TRUE, decreasing=TRUE);
  data.frame(from=v.from[s$ix], to=v.to[s$ix], score=v.score[s$ix], test.value=v.test.value[s$ix], delay=v.delay[s$ix])
}

remove.dup.links <- function(g) {
  # g is the result as returned by infer_grn1
  # For a gene pair x -> y (with direction), if there are more than one entry in g (with different delays),
  #   retain only the one with the highest score (if tie, then larger abs(test.value). If tie again, then smaller delay)
  # Use a naive method
  nr <- max(g$from);
  nc <- max(g$to);
  r <- matrix(rep(0,nr*nc),nrow=nr, ncol=nc); # the (i,j)th entry is the index of the current entry in g to keep for i -> j
  to.keep <- rep(FALSE,nrow(g));
  for(i in 1:nrow(g)) {
    x <- g$from[i];
    y <- g$to[i];
    z <- r[x,y];
    if(z == 0) { # no entry yet, add it
      r[x,y] <- i;
      to.keep[i] <- TRUE;
    } else { # compare
      if((g$score[i] > g$score[z]) ||
         ((g$score[i] == g$score[z]) &&
          ((abs(g$test.value[i]) > abs(g$test.value[z])) || (g$delay[i] < g$delay[z])))) { # better or tie but smaller delay, replace
        r[x,y] <- i;
        to.keep[z] <- FALSE;
        to.keep[i] <- TRUE;
      }
    }
  }
  # done. now take the subset
  subset(g,to.keep)
}
#####

source("pcor.R");

###############################
## Wrapper of the correlation and partial correlation test
test.corr <- function(x,y) {
  # correlation between x and y, in terms of scores
  r <- cor.test(x,y);
  list(score=-log10(r$p.value), test.value=r$estimate)
}

test.p.corr <- function(x,y,z) {
  # partial correlation of x and y, given z (a vector or matrix), in terms of scores
  r <- pcor.test(x,y,z);
  list(score=-log10(r$p.value), test.value=r$estimate)
}
###############################

###############################
# Gaussian kernel probability density to estimate Mutual Information and Conditional Mutual Information
# Formula from Zhang et. al. 2011. "Inferring gene regulatory networks from gene expression data by PC-algorithm based on conditional mutual information". Bioinformatics, 2011.
mutual_info <- function(x,y) {
  # the unconditional one
  # x and y are vectors of the same length
  # Formula: I(X,Y)=(1/2)log(|C(X)|.|C(Y)|/|C(X,Y)|, where C(..) is the covariance matrix, |..| is determinant
  # Here we assume x and y are only vectors, so covariance is easy
  vx <- var(x);
  vy <- var(y);
  cxy <- cov(x,y);
  # since here is the special case, and C(X,Y) is 2 by 2, we use the formula for its determinant directly
  r <- 0.5*log(vx*vy/(vx*vy-cxy*cxy));
  list(score=r, test.value=r)
}

c_mutual_info <- function(x,y,z) {
  # the conditional one
  # x and y should be vectors
  # z can be vector or matrix. In case of matrix, the each column is one vector for one variable
  # Formula: (1/2)log((|C(X,Z)|.|C(Y,Z)|)/(|C(Z)|.|C(X,Y,Z)|))
  cxz <- det(cov(cbind(x,z)));
  cyz <- det(cov(cbind(y,z)));
  cz <- if(is.vector(z)) var(z) else det(cov(z));
  cxyz <- det(cov(cbind(x,y,z)));
  r <- 0.5*log((cxz*cyz)/(cz*cxyz));
  list(score=r, test.value=r)
}
###############################

delayed_exp <- function(r,x,y,d, cs,zs,ts, row.ix=NULL) {
  # r is the expression, row is time points, column is gene
  # x,y are the index of the gene, tentatively x -> y
  # d is the delay of x -> y
  # cs is the vector of indices of the chosen genes in zs
  # zs is the vector of indices of the neighbors to consider, excluded x and y
  # ts is the vector of time delays (relative to x), corresponding to zs
  # row.ix is the vector of row indices of the rows that we really operate on,
  #  i.e. can pick a subset at desired order. If NULL, simply 1:nrow(r)
  # Returns the expressions of x,y and z (the one to be conditioned on) with proper delayed sections
  # return a list: xr, yr and zr, where zr is a matrix
  if(is.null(row.ix)) {row.ix <- 1:nrow(r)}
  ds <- c(0,d,ts[cs]); # delays of: x, y, z
  ds <- ds - min(ds); # to make the delays all non-negative
  L <- length(row.ix) - max(ds); # the overlapping lengths of the delayed expressions
  zsr <- matrix(data=rep(0,length(cs)*L),nrow=L,ncol=length(cs));
  for(i in 1:length(cs)) {
    j <- cs[i];
    dt <- ds[2+i];
    zsr[,i] <- r[row.ix[(1+dt):(dt+L)],zs[j]]; # take the proper delayed section
  }
  # for x
  dt <- ds[1];
  xsr <- r[row.ix[(1+dt):(dt+L)],x];
  # for y
  dt <- ds[2];
  ysr <- r[row.ix[(1+dt):(dt+L)],y];
  # done
  list(xr = xsr, yr = ysr, zr = zsr)
}

gen_combs <- function(lens,cs, ix) {
  # cs is a list indexed by integer, same length as lens
  # each item in cs is a vector of items.
  # ix is a vector of indices into the elements of cs
  # ix[i] means to take an arbitary item from the ith vector in cs
  # Return a matrix, where the columns are the possible combinations of such items
  nc <- prod(lens[ix]); # total number of columns
  nr <- length(ix);
  m <- matrix(rep(cs[[1]][1],nc*nr), nrow=nr, ncol=nc);
  # fill it in by rows, where the first row changes most rapidly, the second less rapidly
  z <- 1; # repeat for how many items, before moving to the next
  for(i in 1:nr) {
    ls <- cs[[ix[i]]];
    L <- lens[ix[i]];
    for(j in 1:nc) {
      # j %/% z gives the item (with possible cycling),
      # then 1+ ((j%/%z)%%lens[i]) gives the proper item
      m[i,j] <- ls[1+(((j-1) %/% z) %% L)];
    }
    z <- z*L;
  }
  m
}

test.c.link <- function(r,st,x,y,n,d, zs,ts, u.neis, test.f,segs=NULL) {
  # r is the expression, row is time points, column is gene
  # st is the score threshold, below which the link is not significant
  # test.f is a function accepting two vectors (of the same length), and
  #  returns a list of (score,test.value), where the score is higher the better,
  #  test.value is whatever value representing the strenth of the test (maybe test statistic, or p-value).
  # segs is a vector of the length of the different segments in r, test.f is applied to 
  #  different segments pieced together.
  # x,y are the index of the gene, tentatively x -> y
  # n is the size of the subset of neighbors to condition on
  # d is the delay of x -> y
  # zs is the vector of indices of the neighbors to consider, excluded x and y
  # ts is the vector of time delays (relative to x), corresponding to zs
  # u.neis is the vector of unique neighbors (the variable id) being conditioned on
  # If the link x -> y is spurious (i.e. score not >= st, i.e. not significant
  #   condition on some subset of variables in u.neis, whose delays come from the possible
  #   multiple links of u.neis), then returns -1, otherwise returns the original delay.
  if(is.null(segs)) {segs <- c(nrow(r))}
  ng <- ncol(r);
  lens <- rep(0,length(u.neis)); # lengths of the list of indices into zs of each neighbor variable
  ixs <- list();
  for(i in 1:length(u.neis)) {
    ixs[[i]] <- (1:length(zs))[zs==u.neis[i]]; # list of indices for one neighbor variable
    lens[i] <- length(ixs[[i]]);
  }
  #
  ns <- combn(1:length(u.neis),n); # each column is a subset of n neighbors
  for(i in 1:ncol(ns)) {
    # generate the combinations in zs for one combination of n neighbors
    cs <- gen_combs(lens,ixs,ns[,i]);
    #
    for(j in 1:ncol(cs)) {
      # first pick the segments
      tmp.x <- NULL;
      tmp.y <- NULL;
      tmp.z <- NULL;
      tmp.L <- 0;
      for(k in 1:length(segs)) {
        w <- delayed_exp(r,x,y,d, cs[,j],zs,ts, (tmp.L+1):(tmp.L+segs[k])); # picking each segment and the conditioned parts
        if(is.null(tmp.x)) {
          tmp.x <- w$xr;
          tmp.y <- w$yr;
          tmp.z <- w$zr;
        } else { # append to the end. Though not too efficient, but simple
          tmp.x <- c(tmp.x, w$xr);
          tmp.y <- c(tmp.y, w$yr);
          tmp.z <- rbind(tmp.z, w$zr);
        }
        tmp.L <- tmp.L + segs[k];
      }
      # take calculate the score
      z <- test.f(tmp.x, tmp.y, tmp.z);
      #
      if(is.numeric(z$score) && !is.na(z$score) && (z$score < st)) {
        return(-1)
      }
    }
  }
  d
}

filter_pcor <- function(r,st,grn, max_n, test.f, pruning=c("all","common"), segs=NULL) {
  # r is as the parameters in infer_grn1
  # st is the score threshold, below which the link is not significant
  # test.f is a function accepting two vectors (of the same length), and
  #  returns a list of (score,test.value), where the score is higher the better,
  #  test.value is whatever value representing the strenth of the test (maybe test statistic, or p-value).
  # grn is a data frame as returned by infer_grn1
  # max_n is the maximum number of neighbors begin conditioned on in partial correlation test
  # If pruning is "all", means for a link, use all neighbors of one vertex to condition on.
  # If pruning is "common", means for a link, condition on the common neighbors of the two vertices.
  # Tries to filter out spurious connections by using partial correlation
  # returns a new data frame by eliminating some edges from grn
  if(nrow(grn) <= 0) {return(grn);}
  pruning <- match.arg(pruning);
  is.common.neis <- pruning == "common";
  #
  ng <- nrow(r);
  n <- 1;
  is_end <- FALSE;
  is_used <- rep(TRUE,nrow(grn)); # to indicate which edges are to be kept
  n_neis <- rep(-1,nrow(grn)); # the number of common neighbors for an edge, -1 means uninitialized
  while((!is_end) && (n <= max_n)) {
    # for each (x,y), check the partial correlation condition on subset of size n
    # common neighbors of one of them
    is_end <- TRUE;
    for(i in nrow(grn):1) {
      if(is_used[i] && (n_neis[i]>=n || n_neis[i]<0)) {
        x <- grn$from[i];
        y <- grn$to[i];
        d <- grn$delay[i];
        # determine the parents and children of x and y
        not.from.x.y <- (grn$from != x) & (grn$from != y) & is_used;
        not.to.x.y <- (grn$to != x) & (grn$to != y) & is_used;
        in.nx <- (grn$to==x) & not.from.x.y;
        out.nx <- (grn$from==x) & not.to.x.y;
        in.nx[i] <- FALSE;
        out.nx[i] <- FALSE;
        # in.nx and out.nx contain the neighbors of x, except itself and y
        in.ny <- (grn$to==y) & not.from.x.y;
        out.ny <- (grn$from==y) & not.to.x.y;
        in.ny[i] <- FALSE;
        out.ny[i] <- FALSE;
        # ny contains the neighbors of y, except itself and x
        # get the distinct neighbor variables for x and y
        neis.x <- c(grn$from[in.nx], grn$to[out.nx]);
        neis.y <- c(grn$from[in.ny], grn$to[out.ny]);
        u.neis.x <- unique(neis.x);
        u.neis.y <- unique(neis.y);
        #
        if(is.common.neis) {
          # calculate the number of common neighbor variables
          u.neis <- intersect(u.neis.x, u.neis.y);
          nxyc <- length(u.neis);
        } else {
          nxc <- length(u.neis.x);
          nyc <- length(u.neis.y);
          nxyc <- min(nxc,nyc);
        }
        n_neis[i] <- nxyc;
        if(nxyc >= n) { # enough neighbors, need to check
          is_end <- FALSE;
          if(is.common.neis) {
            # get the common neighbors with delays relative to x
            zss <- c(neis.x,neis.y);
            tss <- c(-(grn$delay[in.nx]),grn$delay[out.nx], d -(grn$delay[in.ny]), d + grn$delay[out.ny]);
            ixx <- rep(FALSE,length(zss)); # to mark which neighbor is included in the conditioning
            for(j in 1:length(u.neis)) {
              ixx[zss == u.neis[j]] <- TRUE;
            }
            zs <- zss[ixx]; # names of the neighbors to condition on, may have duplicates
            ts <- tss[ixx]; # delays of the neighbors to condition on, relative to x
          } else { # get all the neighbors of one of the vertex
            if(nxc < nyc) {
              # use the neighbors of x
              zs <- neis.x;
              ts <- c(-(grn$delay[in.nx]),grn$delay[out.nx]);
              u.neis <- u.neis.x;
            } else {
              # use the neighbors of y
              zs <- neis.y;
              ts <- d + c(-(grn$delay[in.ny]),grn$delay[out.ny]); # delays relative to x
              u.neis <- u.neis.y;
            }
          }
          # ready
          #cat("Considering ",x," -> ",y,"\td: ",d,"\tscore: ",grn$score[i],"neis: ",zs,"\n", sep=" ");
          d <- test.c.link(r,st,x,y,n,d, zs,ts, u.neis, test.f, segs);
          if(d < 0) {is_used[i] <- FALSE;}
        }
      }
    }
    #
    n <- n+1;
  }
  # done, returns the sub data frame
  subset(grn,is_used)
}

infer_grn <- function(r,st1,st2,lags,max_n, method=c("pcor","mi","other"), u.test.f=NULL,c.test.f=NULL, pruning=c("all","common"), segs=NULL, one.delay=FALSE, no.dup=FALSE) {
  # r, lags parameters same as in infer_grn1
  # st1 is the score threshold for initial links, below which the link is not significant
  # st2 is the score threshold for eliminating links by conditioning, below which the link is not significant
  # method specifies whether to use parital correlation or mutual information for testing significant links,
  #  if the method is "other", then u.test.f and c.test.f are used instead.
  # u.test.f is a function accepting two vectors (of the same length) for testing unconditional significatn links, and
  #  returns a list of (score,test.value), where the score is higher the better,
  #  test.value is whatever value representing the strenth of the test (maybe test statistic, or p-value).
  # c.test.f is a function accepting two vectors, and a vector or matrix,
  #  for testing significant link of the two vectors conditioning on the third.
  # max_n is the maximum number of neighbors begin conditioned on in partial correlation test
  # If pruning is "all", means for a link, use all neighbors of one vertex to condition on.
  # If pruning is "common", means for a link, condition on the common neighbors of the two vertices.
  # segs is a vector of the length of the different segments in r, test.f is applied to each
  #  segment (with each delay), and the scores and test.values are averaged over the segment to be the final value
  # If no.dup is TRUE, then for the final results, for each directed gene pair (x -> y), only the one with
  #  the highest score (if tie, the one with smaller delay) remains.
  # If one.delay is TRUE, then duplicates are removed (as above) after the first stage.
  if(is.null(segs)) {segs <- c(nrow(r))}
  method <- match.arg(method);
  pruning <- match.arg(pruning);
  if(method == "pcor") {u.test.f <- test.corr; c.test.f <- test.p.corr;}
  if(method == "mi") {u.test.f <- mutual_info; c.test.f <- c_mutual_info;}
  z <- infer_grn1(r,st1,lags, u.test.f, segs);
  if(one.delay) {z <- remove.dup.links(z);}
  z <- filter_pcor(r,st2,z, max_n, c.test.f, pruning, segs);
  if(no.dup) remove.dup.links(z) else z
}
