# get_dists() is a function for calling a number of competing methods
# for fitting a hyperplane.  no_hbreg() is for datasets where hbreg
# or other methods from Olive (2017) fail - it runs the other 
# methods.  For some of the more time consuming methods and for
# general hyperplanes, there are limits based on the number of
# dimensions n hard coded in the functions based on observations 
# from experiments.    
# Purpose: to apply competing methods of fitting a hyperplane to 
# data that are contaminated by outliers. 
# Methods:
#  - lm: linear regression
#  - hbreg: composite method from Olive (2017) implemented in 
#          in mpack.txt.  Starts with OLS on the whole dataset.
#          Choose a > 1.  If a*arob is smaller than OLS error, 
#          choose arob; otherwise, if a*bb is smaller than both, 
#          choose bb. a=1.4 is the default.
#  - arob: median ball algorithm to identify outliers before 
#           fitting.
#  - bb: high-breakdown estimator with 10 concentration steps.  Sets
#         are chosen to ensure high breakdown
#  - lts: least trimmed squares regression 
#  - lqs: least quartile of squares regression
#  - alg3: Algorithm 3 from Bertsimas and Mazumder (2014), designed 
#           as a warmup to MIO
#
# Functions, inputs and outputs:
# - fit_lm: fit a linear regression model with the j^th column as the
#           response
#     - inputs: 
#        - j: response column
#        - X: data frame
#        - m: number of normal/non-outlier points
#        - n: number of variables/dimensions
#        - q: percentile for least quantile of squares regression
#     - outputs:
#        - lm_dist: sum of squared orthogonal residuals
#        - lm_norm: normal vector to fitted hyperplane
#        - lm_intercept: intercept for fitted hyperplane
#        - lm_rq: the q^th squared residual, measured orthogonally
#        - lm_lts: the sum of squares of the m smallest residuals, measured orthogonally
# - fit_hbreg: fit an hbreg model, an arob model, and a bb model as 
#              implemented in mpack.txt from Olive (2018).
#    - inputs:
#        - j: response column
#        - X: data frame
#        - m: number of normal/non-outlier points
#        - n: number of variables/dimensions
#        - q: percentile for least quantile of squares regression
#     - outputs:
#        - hbreg_dist: sum of squared orthogonal residuals
#        - hbreg_norm: normal vector to fitted hyperplane
#        - hbreg_intercept: intercept for fitted hyperplane
#        - hbreg_rq: the q^th squared residual
#        - hbreg_lts: the sum of squares of the m smallest residuals, measured orthogonally
#        - arob_dist: sum of squared orthogonal residuals
#        - arob_norm: normal vector to fitted hyperplane
#        - arob_intercept: intercept for fitted hyperplane
#        - arob_rq: the q^th squared residual
#        - arob_lts: the sum of squares of the m smallest residuals, measured orthogonally
#        - bb_dist: sum of squared orthogonal residuals
#        - bb_norm: normal vector to fitted hyperplane
#        - bb_intercept: intercept for fitted hyperplane
#        - bb_rq: the q^th squared residual
#        - bb_lts: the sum of squares of the m smallest residuals, measured orthogonally
# - fit_lts: fit a least trimmed squares regression model with the 
#            j^th column as the response using the method in R.
#     - inputs: 
#        - j: response column
#        - X: data frame
#        - m: number of normal/non-outlier points
#        - n: number of variables/dimensions
#        - q: percentile for least quantile of squares regression
#     - outputs:
#        - lts_dist: sum of squared orthogonal residuals
#        - lts_norm: normal vector to fitted hyperplane
#        - lts_intercept: intercept for fitted hyperplane
#        - lts_rq: the q^th squared residual
#        - lts_lts: the sum of the m smallest residuals, measured orthogonally
# - fit_lqs: fit a least quantile squares regression model with the 
#            j^th column as the response using the method in R.
#     - inputs: 
#        - j: response column
#        - X: data frame
#        - m: number of normal/non-outlier points
#        - n: number of variables/dimensions
#        - q: percentile for least quantile of squares regression
#     - outputs:
#        - lqs_dist: sum of squared orthogonal residuals
#        - lqs_norm: normal vector to fitted hyperplane
#        - lqs_intercept: intercept for fitted hyperplane
#        - lqs_rq: the q^th squared residual
#        - lqs_lts: the sum of squares of the m smallest residuals, measured orthogonally
# - run_alg3: fit a model using Algorithm 3 from Bertsimas and 
#             Mazumder (2014).
#     - inputs:
#        - dataloc: folder containing the data file
#        - srcloc: folder containing implementation of Algorithm 3
#        - fname: data file name
#        - q: percentile for least quantile of squares regression 
#        - dep_var: whether there is a dependent variable or not
#        - resloc: folder for results
#        - mosekloc: location of MOSEK binary 
#        - init_method: method for initializing Algorithm 3 when no
#          response variable is specified.  Options are "PCA" for 
#          principal component analysis and "LP" for elastic LP
#     - outputs: see run_alg3.m
# - get_dists: run lm, hbreg, arob, bb, lts, lqs, RBM, RBM-MIO3, 
#              and algorithm 3 on a dataset
# 
#     - inputs: 
#        - dataloc: folder containing the data file
#        - srcloc: folder containing implementation of Algorithm 3
#        - fname: data file name
#        - q: percentile for least quantile of squares regression 
#        - dep_var: whether there is a dependent variable or not
#        - resloc: folder for results
#        - mosekloc: location of MOSEK binary 
#     - outputs: results are written to various files in resloc 
#                folder
# - no_hbreg: run lm, lts, lqs, RBM, RBM-MIO3,and algorithm 3 on a 
#             dataset
#     - inputs: 
#        - dataloc: folder containing the data file
#        - srcloc: folder containing implementation of Algorithm 3
#        - fname: data file name
#        - q: percentile for least quantile of squares regression 
#        - dep_var: whether there is a dependent variable or not
#        - resloc: folder for results
#        - mosekloc: location of MOSEK binary 
#     - outputs: results are written to various files in resloc 
#                folder
# - hyper_dist_sq: calculate squared residual for a point and a 
#                  hyperplane, measuring along the orthogonal 
#                  projection
#     - inputs:
#        - w: normal vector to hyperplane
#        - b: intercept for hyperplane.  w^T x + b = 0
#        - x: point to be projected on hyperplane  
#     - output: squared distance of point to hyperplane along the 
#               orthogonal projection

## -----------------------------------------------------------------------------------------------------


# fit regression with jth column as the response
fit_lm <- function(j, X, m, n, q) { # for linear regression - lm
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # set j^th column as response, all others are predictors
  
  my_lm <- lm(my_formula, data=X) # fit model
  lm_norm <- my_lm$coefficients[2:length(my_lm$coefficients)] # get normal vector to hyperplane 
  lm_norm[is.na(lm_norm)] <- 0.0
  lm_norm <- append(lm_norm, -1.0, after=j-1) # coefficient for the dependent variable is -1
  lm_intercept <- my_lm$coefficients[1] # get the intercept
  lm_dist <- sum(apply(X[1:m,], 
                           1,  
                           hyper_dist_sq,  
                           w=lm_norm,
                           b=lm_intercept)) # sum of squared orthogonal distances
  lm_res <- apply(X,1,hyper_dist_sq,w=lm_norm,b=lm_intercept) # calculate all squared distances
  lm_rq <- quantile(lm_res, q) # get the q^th largest squared distance
  lm_lts <- sum(sort(lm_res)[1:m]) # sum of m smallest squared residuals, measured orthogonally
  return(list(
    lm_dist=lm_dist, 
    lm_norm=lm_norm,
    lm_intercept=lm_intercept,
    lm_rq=lm_rq,
    lm_lts=lm_lts))
}

fit_hbreg <- function(j, X, m, n, q) { # for hbreg from Olive (2018) mpack.txt
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # set j^th column as response, all others are predictors

  my_hbreg <- hbreg(x=X[,-c(j)], y=X[,j]) # fit model
  hbreg_norm <- my_hbreg$coef[2:length(my_hbreg$coef)] # get normal vector to hyperplane
  hbreg_norm[is.na(hbreg_norm)] <- 0.0
  hbreg_norm <- append(hbreg_norm, -1.0, after=j-1)# coefficient for the dependent variable is -1
  hbreg_intercept <- my_hbreg$coef[1] # get the intercept
  hbreg_dist <- sum(apply(X[1:m,],
                               1,
                               hyper_dist_sq,
                               w=hbreg_norm,
                               b=hbreg_intercept)) # sum of squared orthogonal distances
  hbreg_res <- apply(X,1,hyper_dist_sq,w=hbreg_norm,b=hbreg_intercept) # calculate all squared distances
  hbreg_rq <- quantile(hbreg_res, q) # get the q^th largest squared distance
  hbreg_lts <- sum(sort(hbreg_res)[1:m])# sum of m smallest squared residuals, measured orthogonally

  arob_norm <- my_hbreg$arobcoef[2:length(my_hbreg$arobcoef)]# get normal vector to hyperplane
  arob_norm[is.na(arob_norm)] <- 0.0
  arob_norm <- append(arob_norm, -1.0, after=j-1)# coefficient for the dependent variable is -1
  arob_intercept <- my_hbreg$arobcoef[1]# get the intercept
  arob_dist <- sum(apply(X[1:m,],
                               1,
                               hyper_dist_sq,
                               w=arob_norm,
                               b=arob_intercept))# sum of squared orthogonal distances
  arob_res <- apply(X,1,hyper_dist_sq,w=arob_norm,b=arob_intercept)# calculate all squared distances
  arob_rq <- quantile(arob_res, q)# get the q^th largest squared distance
  arob_lts <- sum(sort(arob_res)[1:m])# sum of m smallest squared residuals, measured orthogonally
  bb_norm <- my_hbreg$bbcoef[2:length(my_hbreg$bbcoef)]# get normal vector to hyperplane
  bb_norm[is.na(bb_norm)] <- 0.0
  bb_norm <- append(bb_norm, -1.0, after=j-1)# coefficient for the dependent variable is -1
  bb_intercept <- my_hbreg$bbcoef[1]# get the intercept
  bb_dist <- sum(apply(X[1:m,],
                               1,
                               hyper_dist_sq,
                               w=bb_norm,
                               b=bb_intercept)) # sum of squared orthogonal distances

  bb_res <- apply(X,1,hyper_dist_sq,w=bb_norm,b=bb_intercept)# calculate all squared distances
  bb_rq <- quantile(bb_res, q)# get the q^th largest squared distance
  bb_lts <- sum(sort(bb_res)[1:m])# sum of m smallest squared residuals, measured orthogonally

  return(list(
              hbreg_dist=hbreg_dist,
              hbreg_norm=hbreg_norm,
              hbreg_intercept=hbreg_intercept,
              hbreg_rq=hbreg_rq,
              hbreg_lts=hbreg_lts,
              arob_dist=arob_dist,
              arob_norm=arob_norm,
              arob_intercept=arob_intercept,
              arob_rq=arob_rq,
              arob_lts=arob_lts,
              bb_dist=bb_dist,
              bb_norm=bb_norm,
              bb_intercept=bb_intercept,
              bb_rq=bb_rq,
              bb_lts=bb_lts
            )
        )
}

# fit regression with jth column as the response
fit_lts <- function(j, X, m, n, q) {
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep=""))
  my_lts <- lqs(my_formula, data=X, method="lts", quantile=floor(q*m)) 
  lts_norm <- my_lts$coefficients[2:length(my_lts$coefficients)]
  lts_norm[is.na(lts_norm)] <- 0.0
  lts_norm <- append(lts_norm, -1.0, after=j-1)
  lts_intercept <- my_lts$coefficients[1]
  lts_dist <- sum(apply(X[1:m,],
                       1,
                       hyper_dist_sq,
                       w=lts_norm,
                       b=lts_intercept))
  lts_res <- apply(X,1,hyper_dist_sq,w=lts_norm,b=lts_intercept)
  lts_rq <- quantile(lts_res, q)
  lts_lts <- sum(sort(lts_res)[1:m])

  return(list(
              lts_dist=lts_dist,
              lts_norm=lts_norm,
              lts_intercept=lts_intercept,
              lts_rq=lts_rq,
              lts_lts=lts_lts
            )
        )
}

# fit regression with jth column as the response
fit_lqs <- function(j, X, m, n, q) {
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep=""))

  my_lqs <- lqs(my_formula, data=X, method="lqs", quantile=floor(q*m)) 
  lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)]
  lqs_norm[is.na(lqs_norm)] <- 0.0
  lqs_norm <- append(lqs_norm, -1.0, after=j-1)
  lqs_intercept <- my_lqs$coefficients[1]
  lqs_dist <- sum(apply(X[1:m,],
                       1,
                       hyper_dist_sq,
                       w=lqs_norm,
                       b=lqs_intercept))
  
  lqs_res <- apply(X,1,hyper_dist_sq,w=lqs_norm,b=lqs_intercept)
  lqs_rq <- quantile(lqs_res, q)
  lqs_lts <- sum(sort(lqs_res)[1:m])

  return(list(
              lqs_dist=lqs_dist,
              lqs_norm=lqs_norm,
              lqs_intercept=lqs_intercept,
              lqs_rq=lqs_rq,
              lqs_lts=lqs_lts
            )
        )
}


run_alg3 <- function(dataloc, srcloc, fname, q, dep_var, resloc, mosekloc, init_method) {
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname)
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")
  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE)
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")
  q <- floor(q*nrow(X))
  print(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','",init_method,"')", sep=""))
  run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','",init_method,"')", sep=""))
}

## -----------------------------------------------------------------------------------------------------
get_dists <- function(dataloc, srcloc, fname, q, dep_var, timelimit, resloc, mosekloc) {
  #loc <- 'vary_everything'
  #fname <- 'm102n47m_outliers5num_clust1same_sideFALSEoutlier_dist1000i848.csv'
  #fname <- "m650n10m_outliers350i0.csv"
  #loc <- 'rd'
  #library(MASS)
  #source("../2021_02_hyper_data/olive/mpack.txt")
  #q <- 0.64
  print("srcloc")
  print(srcloc)
  
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname)
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")
  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE)
  if (dep_var == TRUE) {
      lm_time <- proc.time()
      my_lm <- lm(V1 ~ ., data=X)  
      lm_dist <- sum(my_lm$residuals[1:m]^2)
      lm_rq <- quantile(abs(my_lm$residuals), q)
      lm_lts <- sum(sort(my_lm$residuals^2)[1:m])
      lm_time <- proc.time() - lm_time

      hbreg_time <- proc.time()
      my_hbreg <- hbreg(x=X[,-1], y=X[,1])
      hbreg_norm <- my_hbreg$coef[2:length(my_hbreg$coef)]
      hbreg_res <- X[,1] - as.matrix(X[,2:n]) %*% hbreg_norm - my_hbreg$coef[1]
      hbreg_dist <- sum(hbreg_res[1:m]^2)
      hbreg_rq <- quantile(abs(hbreg_res), q)
      hbreg_lts <- sum(sort(hbreg_res^2)[1:m])
      hbreg_time <- proc.time() - hbreg_time

      arob_time <- proc.time()
      arob_norm <- my_hbreg$arobcoef[2:length(my_hbreg$arobcoef)]
      arob_res <- X[,1] - as.matrix(X[,2:n]) %*% arob_norm - my_hbreg$arobcoef[1]
      arob_dist <- sum(arob_res[1:m]^2)
      arob_rq <- quantile(abs(arob_res), q)
      arob_lts <- sum(sort(arob_res^2)[1:m])
      arob_time <- proc.time() - arob_time

      bb_time <- proc.time() 
      bb_norm <- my_hbreg$bbcoef[2:length(my_hbreg$bbcoef)]
      bb_res <- X[,1] - as.matrix(X[,2:n]) %*% bb_norm - my_hbreg$bbcoef[1]
      bb_dist <- sum(bb_res[1:m]^2)
      bb_rq <- quantile(abs(bb_res), q)
      bb_lts <- sum(sort(bb_res^2)[1:m])
      bb_time <- proc.time() - bb_time

      lts_time <- proc.time() 
      my_lts <- lqs(V1 ~ ., data=X, method="lts", quantile=floor(q*m)) 
      lts_norm <- my_lts$coefficients[2:length(my_lts$coefficients)]
      lts_res <- X[,1] - as.matrix(X[,2:n]) %*% lts_norm - my_lts$coefficients[1]
      lts_dist <- sum(lts_res[1:m]^2)
      lts_rq <- quantile(abs(lts_res), q)
      lts_lts <- sum(sort(lts_res^2)[1:m])
      lts_time <- proc.time() - lts_time

      lqs_time <- proc.time() 
      my_lqs <- lqs(V1 ~ ., data=X, method="lqs", quantile=floor(q*m)) 
      lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)]
      lqs_res <- X[,1] - as.matrix(X[,2:n]) %*% lqs_norm - my_lqs$coefficients[1]
      lqs_dist <- sum(lqs_res[1:m]^2)
      lqs_rq <- quantile(abs(lqs_res), q)
      lqs_beta <- matrix(c(-1.0, my_lqs$coefficients),
                          nrow=length(lqs_norm)+2,
                          ncol=1) # beta_1 = -1 - the response
      lqs_lts <- sum(sort(lqs_res^2)[1:m])
      lqs_time <- proc.time() - lqs_time

      lm_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lm", lm_dist, lm_rq, lm_time[1], lm_time[2], lm_time[3], lm_lts, sep=",")
      write(lm_out, file=paste(resloc, "/lmi", i, ".csv", sep=""), ncol=length(lm_out))
      hbreg_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "hbreg", hbreg_dist, hbreg_rq, hbreg_time[1], hbreg_time[2], hbreg_time[3], hbreg_lts, sep=",")
      write(hbreg_out, file=paste(resloc, "/hbregi", i, ".csv", sep=""), ncol=length(hbreg_out))
      arob_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "arob", arob_dist, arob_rq, arob_time[1], arob_time[2], arob_time[3], arob_lts, sep=",")
      write(arob_out, file=paste(resloc, "/arobi", i, ".csv", sep=""), ncol=length(arob_out))
      bb_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "bb", bb_dist, bb_rq, bb_time[1], bb_time[2], bb_time[3], bb_lts, sep=",")
      write(bb_out, file=paste(resloc, "/bbi", i, ".csv", sep=""), ncol=length(bb_out))
      lts_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lts", lts_dist, lts_rq, lts_time[1], lts_time[2], lts_time[3], lts_lts, sep=",")
      write(lts_out, file=paste(resloc, "/ltsi", i, ".csv", sep=""), ncol=length(lts_out))
      lqs_out <- paste(paste(resloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lqs", lqs_dist, lqs_rq, lqs_time[1], lqs_time[2], lqs_time[3], lqs_lts, sep=",")
      write(lqs_out, file=paste(resloc, "/lqsi", i, ".csv", sep=""), ncol=length(lqs_out))
 
  } else {
      lm_time <- proc.time()
      lm_results <- lapply(1:n, fit_lm, X, m, n, q)
      lm_dists <- sapply(lm_results, function(x) x$lm_dist)
      lm_dist <- min(lm_dists)
      lm_norms <- sapply(lm_results, function(x) x$lm_norm)
      lm_intercepts <- sapply(lm_results, function(x) x$lm_intercept)
      lm_norm <- lm_norms[,which(lm_dists == min(lm_dists))]
      lm_intercept <- lm_intercepts[which(lm_dists == min(lm_dists))]
      lm_rqs <- sapply(lm_results, function(x) x$lm_rq)
      lm_rq <- lm_rqs[which(lm_dists == min(lm_dists))]
      lm_ltss <- sapply(lm_results, function(x) x$lm_lts)
      lm_lts <- lm_ltss[which(lm_dists == min(lm_dists))]
      lm_time <- proc.time() - lm_time

      lm_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lm", lm_dist, lm_rq, lm_time[1], lm_time[2], lm_time[3], lm_lts, sep=",")
      print(lm_out)
      write(lm_out, file=paste(resloc, "/lmi", i, ".csv", sep=""))

      
      if (n <= 350) {
          hbreg_time <- proc.time()
          hbreg_results <- lapply(1:n, fit_hbreg, X, m, n, q)
          hbreg_dists <- sapply(hbreg_results, function(x) x$hbreg_dist)
          hbreg_dist <- min(hbreg_dists)
          hbreg_norms <- sapply(hbreg_results, function(x) x$hbreg_norm)
          hbreg_intercepts <- sapply(hbreg_results, function(x) x$hbreg_intercept)
          hbreg_norm <- hbreg_norms[,which(hbreg_dists == min(hbreg_dists))]
          hbreg_intercept <- hbreg_intercepts[which(hbreg_dists == min(hbreg_dists))]
          hbreg_rqs <- sapply(hbreg_results, function(x) x$hbreg_rq)
          hbreg_rq <- hbreg_rqs[which(hbreg_dists == min(hbreg_dists))]
          hbreg_ltss <- sapply(hbreg_results, function(x) x$hbreg_lts)
          hbreg_lts <- hbreg_ltss[which(hbreg_dists == min(hbreg_dists))]
          
          arob_dists <- sapply(hbreg_results, function(x) x$arob_dist)
          arob_dist <- min(arob_dists)
          arob_norms <- sapply(hbreg_results, function(x) x$arob_norm)
          arob_intercepts <- sapply(hbreg_results, function(x) x$arob_intercept)
          arob_norm <- arob_norms[,which(arob_dists == min(arob_dists))]
          arob_intercept <- arob_intercepts[which(arob_dists == min(arob_dists))]
          arob_rqs <- sapply(hbreg_results, function(x) x$arob_rq)
          arob_rq <- arob_rqs[which(arob_dists == min(arob_dists))]
          arob_ltss <- sapply(hbreg_results, function(x) x$arob_lts)
          arob_lts <- arob_ltss[which(arob_dists == min(arob_dists))]

          bb_dists <- sapply(hbreg_results, function(x) x$bb_dist)
          bb_norms <- sapply(hbreg_results, function(x) x$bb_norm)
          bb_dist <- min(bb_dists)
          bb_intercepts <- sapply(hbreg_results, function(x) x$bb_intercept)
          bb_norm <- bb_norms[,which(bb_dists == min(bb_dists))]
          bb_intercept <- bb_intercepts[which(bb_dists == min(bb_dists))]
          bb_rqs <- sapply(hbreg_results, function(x) x$bb_rq)
          bb_rq <- bb_rqs[which(bb_dists == min(bb_dists))]
          bb_ltss <- sapply(hbreg_results, function(x) x$bb_lts)
          bb_lts <- bb_ltss[which(bb_dists == min(bb_dists))]

          hbreg_time <- proc.time() - hbreg_time
          hbreg_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "hbreg", hbreg_dist, hbreg_rq, hbreg_time[1], hbreg_time[2], hbreg_time[3], hbreg_lts, sep=",")
          write(hbreg_out, file=paste(resloc, "/hbregi", i, ".csv", sep=""))
          arob_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "arob", arob_dist, arob_rq, hbreg_time[1], hbreg_time[2], hbreg_time[3], arob_lts, sep=",")
          write(arob_out, file=paste(resloc, "/arobi", i, ".csv", sep=""))
          bb_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "bb", bb_dist, bb_rq, hbreg_time[1], hbreg_time[2], hbreg_time[3], bb_lts, sep=",")
          write(bb_out, file=paste(resloc, "/bbi", i, ".csv", sep=""))
      }

      if (n <= 150) { # number of variables
          #write(c(lm_dist, hbreg_dist, arob_dist, bb_dist, lts_dist, lqs_dist), file=paste(loc,"/",fname,"dist", sep=""),ncolumns=6)
          lts_time <- proc.time() 
          lts_results <- lapply(1:n, fit_lts, X, m, n, q)
          lts_dists <- sapply(lts_results, function(x) x$lts_dist)
          lts_dist <- min(lts_dists)
          lts_norms <- sapply(lts_results, function(x) x$lts_norm)
          lts_intercepts <- sapply(lts_results, function(x) x$lts_intercept)
          lts_norm <- lts_norms[,which(lts_dists == min(lts_dists))]
          lts_intercept <- lts_intercepts[which(lts_dists == min(lts_dists))]
          lts_rqs <- sapply(lts_results, function(x) x$lts_rq)
          lts_rq <- lts_rqs[which(lts_dists == min(lts_dists))]
          lts_ltss <- sapply(lts_results, function(x) x$lts_lts)
          lts_lts <- lts_ltss[which(lts_dists == min(lts_dists))]
          lts_time <- proc.time() - lts_time

          lqs_time <- proc.time() 
          lqs_results <- lapply(1:n, fit_lqs, X, m, n, q)
          lqs_dists <- sapply(lqs_results, function(x) x$lqs_dist)
          lqs_dist <- min(lqs_dists)
          lqs_norms <- sapply(lqs_results, function(x) x$lqs_norm)
          lqs_intercepts <- sapply(lqs_results, function(x) x$lqs_intercept)
          lqs_norm <- lqs_norms[,which(lqs_dists == min(lqs_dists))]
          lqs_intercept <- lqs_intercepts[which(lqs_dists == min(lqs_dists))]
          lqs_rqs <- sapply(lqs_results, function(x) x$lqs_rq)
          lqs_rq <- lqs_rqs[which(lqs_dists == min(lqs_dists))]

          lqs_beta <- matrix(c(lqs_intercept, lqs_norm),
                             nrow=length(lqs_norm)+1,
                             ncol=1)
          lqs_beta <- (lqs_beta/lqs_beta[1,1])*n  # normalize so that beta_0 = n
          #lqs_beta <- lqs_beta[,1]
          lqs_ltss <- sapply(lqs_results, function(x) x$lqs_lts)
          lqs_lts <- lqs_ltss[which(lqs_dists == min(lqs_dists))]
          lqs_time <- proc.time() - lqs_time
      
          print(resloc)

          lts_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lts", lts_dist, lts_rq, lts_time[1], lts_time[2], lts_time[3], lts_lts, sep=",")
          write(lts_out, file=paste(resloc, "/ltsi", i, ".csv", sep=""))
          lqs_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lqs", lqs_dist, lqs_rq, lqs_time[1], lqs_time[2], lqs_time[3], lqs_lts, sep=",")
          write(lqs_out, file=paste(resloc, "/lqsi", i, ".csv", sep=""))
      } else {
          lqs_beta <- -10
      }
  }
  print("lqs_beta")
  print(lqs_beta)
  make_lqs_beta <- rmat_to_matlab_mat(lqs_beta, matname="lqs_beta")
  
#  write(c(lm_rq, hbreg_rq, arob_rq, bb_rq, lts_rq, lqs_rq), file=paste(loc,"/",fname,"rq", sep=""),ncolumns=6)

  q <- floor(q*nrow(X))
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")
  print(paste("addpath('", srcloc, "');", sep=""))

  # need to run RBM and get RBM start
  # run rbm when dep_var is FALSE
  if (dep_var == TRUE) {
    run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"','qout',",m, ",true,'", "cbmio3", "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "','PCA')", sep=""))
  } else {
    print(add_path)
    run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"','qout',",m, ",false,'", "cbmio3", "','", resloc, "',", timelimit, ")", sep="")) # dep_var = FALSE
    if (n <= 500) {
      run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','PCA')", sep=""))
    }
  }

  #write(c(lm_norm, lm_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1)
  #write(c(hbreg_norm, hbreg_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1, append=TRUE)
  #write(c(arob_norm, arob_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1, append=TRUE)
  #write(c(bb_norm, bb_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1, append=TRUE)
  #write(c(lts_norm, lts_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1, append=TRUE)
  #write(c(lqs_norm, lqs_intercept), file=paste(loc,"/hp/",fname,"hp", sep=""),ncolumns=n+1, append=TRUE)

}

no_hbreg <- function(dataloc, srcloc, fname, q, dep_var, timelimit, resloc, mosekloc) {
  #loc <- 'vary_everything'
  #fname <- 'm102n47m_outliers5num_clust1same_sideFALSEoutlier_dist1000i848.csv'
  #fname <- "m650n10m_outliers350i0.csv"
  #loc <- 'rd'
  #library(MASS)
  #source("../2021_02_hyper_data/olive/mpack.txt")
  #q <- 0.64
  
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname)
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")
  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE)
  if (dep_var == TRUE) {
      lm_time <- proc.time()
      my_lm <- lm(V1 ~ ., data=X)  
      lm_dist <- sum(my_lm$residuals[1:m]^2)
      lm_rq <- quantile(abs(my_lm$residuals), q)
      lm_lts <- sum(sort(my_lm$residuals^2)[1:m])
      lm_time <- proc.time() - lm_time

      lts_time <- proc.time() 
      my_lts <- lqs(V1 ~ ., data=X, method="lts", quantile=floor(q*m)) 
      lts_norm <- my_lts$coefficients[2:length(my_lts$coefficients)]
      lts_res <- X[,1] - as.matrix(X[,2:n]) %*% lts_norm - my_lts$coefficients[1]
      lts_dist <- sum(lts_res[1:m]^2)
      lts_rq <- quantile(abs(lts_res), q)
      lts_lts <- sum(sort(lts_res^2)[1:m])
      lts_time <- proc.time() - lts_time

      lqs_time <- proc.time() 
      my_lqs <- lqs(V1 ~ ., data=X, method="lqs", quantile=floor(q*m)) 
      lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)]
      lqs_res <- X[,1] - as.matrix(X[,2:n]) %*% lqs_norm - my_lqs$coefficients[1]
      lqs_dist <- sum(lqs_res[1:m]^2)
      lqs_rq <- quantile(abs(lqs_res), q)
      lqs_beta <- matrix(c(-1.0, my_lqs$coefficients),
                          nrow=length(lqs_norm)+2,
                          ncol=1) # beta_1 = -1 - the response
      lqs_lts <- sum(sort(lqs_res^2)[1:m])
      lqs_time <- proc.time() - lqs_time

      lm_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lm", lm_dist, lm_rq, lm_time[1], lm_time[2], lm_time[3], lm_lts, sep=",")
      write(lm_out, file=paste(resloc, "/lmi", i, ".csv", sep=""), ncol=length(lm_out))
      lts_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lts", lts_dist, lts_rq, lts_time[1], lts_time[2], lts_time[3], lts_lts, sep=",")
      write(lts_out, file=paste(resloc, "/ltsi", i, ".csv", sep=""), ncol=length(lts_out))
      lqs_out <- paste(paste(resloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lqs", lqs_dist, lqs_rq, lqs_time[1], lqs_time[2], lqs_time[3], lqs_lts, sep=",")
      write(lqs_out, file=paste(resloc, "/lqsi", i, ".csv", sep=""), ncol=length(lqs_out))
 
  } else {
      lm_time <- proc.time()
      lm_results <- lapply(1:n, fit_lm, X, m, n, q)
      lm_dists <- sapply(lm_results, function(x) x$lm_dist)
      lm_dist <- min(lm_dists)
      lm_norms <- sapply(lm_results, function(x) x$lm_norm)
      lm_intercepts <- sapply(lm_results, function(x) x$lm_intercept)
      lm_norm <- lm_norms[,which(lm_dists == min(lm_dists))]
      lm_intercept <- lm_intercepts[which(lm_dists == min(lm_dists))]
      lm_rqs <- sapply(lm_results, function(x) x$lm_rq)
      lm_rq <- lm_rqs[which(lm_dists == min(lm_dists))]
      lm_ltss <- sapply(lm_results, function(x) x$lm_lts)
      lm_lts <- lm_ltss[which(lm_dists == min(lm_dists))]
      lm_time <- proc.time() - lm_time
      lm_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lm", lm_dist, lm_rq, lm_time[1], lm_time[2], lm_time[3], lm_lts, sep=",")
      print(lm_out)
      write(lm_out, file=paste(resloc, "/lmi", i, ".csv", sep=""))
      
      if (n <= 150) { # number of variables
          #write(c(lm_dist, hbreg_dist, arob_dist, bb_dist, lts_dist, lqs_dist), file=paste(loc,"/",fname,"dist", sep=""),ncolumns=6)
          
          lts_time <- proc.time() 
          lts_results <- lapply(1:n, fit_lts, X, m, n, q)
          lts_dists <- sapply(lts_results, function(x) x$lts_dist)
          lts_dist <- min(lts_dists)
          lts_norms <- sapply(lts_results, function(x) x$lts_norm)
          lts_intercepts <- sapply(lts_results, function(x) x$lts_intercept)
          lts_norm <- lts_norms[,which(lts_dists == min(lts_dists))]
          lts_intercept <- lts_intercepts[which(lts_dists == min(lts_dists))]
          lts_rqs <- sapply(lts_results, function(x) x$lts_rq)
          lts_rq <- lts_rqs[which(lts_dists == min(lts_dists))]
          lts_ltss <- sapply(lts_results, function(x) x$lts_lts)
          lts_lts <- lts_ltss[which(lts_dists == min(lts_dists))]
          lts_time <- proc.time() - lts_time

          lqs_time <- proc.time() 
          lqs_results <- lapply(1:n, fit_lqs, X, m, n, q)
          lqs_dists <- sapply(lqs_results, function(x) x$lqs_dist)
          lqs_dist <- min(lqs_dists)
          lqs_norms <- sapply(lqs_results, function(x) x$lqs_norm)
          lqs_intercepts <- sapply(lqs_results, function(x) x$lqs_intercept)
          lqs_norm <- lqs_norms[,which(lqs_dists == min(lqs_dists))]
          lqs_intercept <- lqs_intercepts[which(lqs_dists == min(lqs_dists))]
          lqs_rqs <- sapply(lqs_results, function(x) x$lqs_rq)
          lqs_rq <- lqs_rqs[which(lqs_dists == min(lqs_dists))]

          lqs_beta <- matrix(c(lqs_intercept, lqs_norm),
                             nrow=length(lqs_norm)+1,
                             ncol=1)
          lqs_beta <- (lqs_beta/lqs_beta[1,1])*n  # normalize so that beta_0 = n
          #lqs_beta <- lqs_beta[,1]
          lqs_ltss <- sapply(lqs_results, function(x) x$lqs_lts)
          lqs_lts <- lqs_ltss[which(lqs_dists == min(lqs_dists))]
          lqs_time <- proc.time() - lqs_time
      
          print(resloc)

          lts_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lts", lts_dist, lts_rq, lts_time[1], lts_time[2], lts_time[3], lts_lts, sep=",")
          write(lts_out, file=paste(resloc, "/ltsi", i, ".csv", sep=""))
          lqs_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lqs", lqs_dist, lqs_rq, lqs_time[1], lqs_time[2], lqs_time[3], lqs_lts, sep=",")
          write(lqs_out, file=paste(resloc, "/lqsi", i, ".csv", sep=""))
      } else {
          lqs_beta <- -10
      }
  }
  print("lqs_beta")
  print(lqs_beta)
  make_lqs_beta <- rmat_to_matlab_mat(lqs_beta, matname="lqs_beta")
  

  q <- floor(q*nrow(X))
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")
  print(paste("addpath('", srcloc, "');", sep=""))

  # need to run RBM and get RBM start
  # run rbm when dep_var is FALSE
  if (dep_var == TRUE) {
    run_matlab_code(paste(add_path, " ", "rbm(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "')", sep=""))
    run_matlab_code(paste(add_path, " ", make_lqs_beta, " ",  "rbmmio3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "',", timelimit, ",lqs_beta)", sep=""))
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "','PCA')", sep=""))
  } else {
    run_matlab_code(paste(add_path, " ", "rbm(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "')", sep=""))
    run_matlab_code(paste(add_path, " ", make_lqs_beta, " ",  "rbmmio3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "',", timelimit, ",lqs_beta)", sep=""))
    if (n <= 500) {
      run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','PCA')", sep=""))
    }
  }

}

run_cb <- function(dataloc, srcloc, fname, q, dep_var, resloc, mosekloc) {
  print("srcloc")
  print(srcloc)
  
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname)
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")

  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE)
  q <- floor(q*nrow(X))
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")
  print(paste("addpath('", srcloc, "');", sep=""))

  # need to run RBM and get RBM start
  # run rbm when dep_var is FALSE
  if (dep_var == TRUE) {
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "')", sep=""))
  } else {
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "')", sep=""))
  }
}


## -----------------------------------------------------------------------------------------------------
hyper_dist_sq <- function(w,b,x) {
  (sum(w*x) + b)^2/sum(w*w)
}


