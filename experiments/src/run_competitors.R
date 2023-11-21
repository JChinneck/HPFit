# get_dists() is a function for calling a number of competing methods
# for fitting a hyperplane.  The goal is to compare methods for 
# minimizing the trimmed sum of errors.  
# no_hbreg() is for datasets where hbreg
# or other methods from Olive (2017) fail - it runs the other 
# methods.  For some of the more time consuming methods and for
# general hyperplanes, there are limits based on the number of
# dimensions n hard coded in the functions based on observations 
# from experiments.    
# Other supporting functions include:
# fit_lm()    fit an OLS model with a specified variable as dependent
# fit_hbreg() fit hbreg, arob, bb with a specified variable as 
#             dependent
# fit_lts()   fit an lts model with a specified variable as dependent
# fit_lqs()   fit an lqs model with a specified variable as dependent
# lqs_gamma() fit an lqs model with a specified variable as dependent
# hyper_dist_sq() get the distance to a hyperplane

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
# - get_dists: run lm, hbreg, arob, bb, lts, lqs, CB, 
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
# - no_hbreg: run lm, lts, lqs, CB, and algorithm 3 on a 
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

  arob_time <- my_hbreg$arob_time
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
 
  bb_time <- my_hbreg$bb_time
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
              arob_time=arob_time,
              arob_dist=arob_dist,
              arob_norm=arob_norm,
              arob_intercept=arob_intercept,
              arob_rq=arob_rq,
              arob_lts=arob_lts,
              bb_time=bb_time,
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

# fit regression with jth column as the response
lqs_gamma <- function(j, X, m, n, q) {
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep=""))

  my_lqs <- lqs(my_formula, data=X, method="lqs", quantile=floor(q*m)) 
  lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)]
  lqs_norm[is.na(lqs_norm)] <- 0.0
  lqs_norm <- append(lqs_norm, -1.0, after=j-1)
  lqs_intercept <- my_lqs$coefficients[1]
  lqs_beta <- matrix(c(lqs_intercept, lqs_norm),
                         nrow=length(lqs_norm)+1,
                         ncol=1)
  lqs_beta <- (lqs_beta/lqs_beta[1,1])*n  # normalize so that beta_0 = n
  
  lqs_res <- apply(X,1,hyper_dist_sq,w=lqs_beta[2:nrow(lqs_beta),1],b=n)
  lqs_rq <- quantile(lqs_res, q)
  lqs_lts <- sum(sort(lqs_res)[1:m])
  lqs_dist <- lqs_rq

  return(list(
              lqs_dist=lqs_dist,
              lqs_norm=lqs_norm,
              lqs_intercept=lqs_intercept,
              lqs_rq=lqs_rq,
              lqs_lts=lqs_lts
            )
        )
}

# fit M regression with jth column as the response
fit_mh <- function(j, X, m, n, q) { # for M regression - rlm
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # set j^th column as response, all others are predictors
  
  my_mh <- rlm(my_formula, data=X, method="M") # fit model
  mh_norm <- my_mh$coefficients[2:length(my_mh$coefficients)] # get normal vector to hyperplane 
  mh_norm[is.na(mh_norm)] <- 0.0
  mh_norm <- append(mh_norm, -1.0, after=j-1) # coefficient for the dependent variable is -1
  mh_intercept <- my_mh$coefficients[1] # get the intercept
  mh_dist <- sum(apply(X[1:m,], 
                           1,  
                           hyper_dist_sq,  
                           w=mh_norm,
                           b=mh_intercept)) # sum of squared orthogonal distances
  mh_res <- apply(X,1,hyper_dist_sq,w=mh_norm,b=mh_intercept) # calculate all squared distances
  mh_rq <- quantile(mh_res, q) # get the q^th largest squared distance
  mh_lts <- sum(sort(mh_res)[1:m]) # sum of m smallest squared residuals, measured orthogonally
  return(list(
    mh_dist=mh_dist, 
    mh_norm=mh_norm,
    mh_intercept=mh_intercept,
    mh_rq=mh_rq,
    mh_lts=mh_lts))
}

# fit MM regression with jth column as the response
fit_mm <- function(j, X, m, n, q) { # for MM regression - lmRob
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # set j^th column as response, all others are predictors
  
  my_mm <- lmRob(my_formula, data=X) # fit model
  mm_norm <- my_mm$coefficients[2:length(my_mm$coefficients)] # get normal vector to hyperplane 
  mm_norm[is.na(mm_norm)] <- 0.0
  mm_norm <- append(mm_norm, -1.0, after=j-1) # coefficient for the dependent variable is -1
  mm_intercept <- my_mm$coefficients[1] # get the intercept
  mm_dist <- sum(apply(X[1:m,], 
                           1,  
                           hyper_dist_sq,  
                           w=mm_norm,
                           b=mm_intercept)) # sum of squared orthogonal distances
  mm_res <- apply(X,1,hyper_dist_sq,w=mm_norm,b=mm_intercept) # calculate all squared distances
  mm_rq <- quantile(mm_res, q) # get the q^th largest squared distance
  mm_lts <- sum(sort(mm_res)[1:m]) # sum of m smallest squared residuals, measured orthogonally
  return(list(
    mm_dist=mm_dist, 
    mm_norm=mm_norm,
    mm_intercept=mm_intercept,
    mm_rq=mm_rq,
    mm_lts=mm_lts))
}

# fit REWLSE regression with jth column as the response
fit_rewlse <- function(j, X, m, n, q) { # for REWLSE regression - lmRob
  cat(j, " ")
  #j<-1
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # set j^th column as response, all others are predictors
  
  my_rewlse <- lmRob(my_formula, data=X, control=lmRob.control(final.alg="Adaptive")) # fit model
  rewlse_norm <- my_rewlse$coefficients[2:length(my_rewlse$coefficients)] # get normal vector to hyperplane 
  rewlse_norm[is.na(rewlse_norm)] <- 0.0
  rewlse_norm <- append(rewlse_norm, -1.0, after=j-1) # coefficient for the dependent variable is -1
  rewlse_intercept <- my_rewlse$coefficients[1] # get the intercept
  rewlse_dist <- sum(apply(X[1:m,], 
                           1,  
                           hyper_dist_sq,  
                           w=rewlse_norm,
                           b=rewlse_intercept)) # sum of squared orthogonal distances
  rewlse_res <- apply(X,1,hyper_dist_sq,w=rewlse_norm,b=rewlse_intercept) # calculate all squared distances
  rewlse_rq <- quantile(rewlse_res, q) # get the q^th largest squared distance
  rewlse_lts <- sum(sort(rewlse_res)[1:m]) # sum of m smallest squared residuals, measured orthogonally
  return(list(
    rewlse_dist=rewlse_dist, 
    rewlse_norm=rewlse_norm,
    rewlse_intercept=rewlse_intercept,
    rewlse_rq=rewlse_rq,
    rewlse_lts=rewlse_lts))
}


run_alg3 <- function(dataloc, srcloc, fname, q, dep_var, resloc, mosekloc, init_method) {
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname)
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")
  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE)

  q <- floor(q*nrow(X))
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")

  if (dep_var == TRUE) {
    dep_var_matlab <- "true" 
  } else{
    dep_var_matlab <- "false"
  }

  print(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m,",",dep_var_matlab,",'", resloc, "','",init_method,"')", sep=""))
  run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m,",",dep_var_matlab,",'", resloc, "','",init_method,"')", sep=""))
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
      lqs_beta <- -10
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

      arob_time <- my_hbreg$arob_time
      arob_norm <- my_hbreg$arobcoef[2:length(my_hbreg$arobcoef)]
      arob_res <- X[,1] - as.matrix(X[,2:n]) %*% arob_norm - my_hbreg$arobcoef[1]
      arob_dist <- sum(arob_res[1:m]^2)
      arob_rq <- quantile(abs(arob_res), q)
      arob_lts <- sum(sort(arob_res^2)[1:m])

      bb_time <- my_hbreg$bb_time
      bb_norm <- my_hbreg$bbcoef[2:length(my_hbreg$bbcoef)]
      bb_res <- X[,1] - as.matrix(X[,2:n]) %*% bb_norm - my_hbreg$bbcoef[1]
      bb_dist <- sum(bb_res[1:m]^2)
      bb_rq <- quantile(abs(bb_res), q)
      bb_lts <- sum(sort(bb_res^2)[1:m])

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
      lqs_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "lqs", lqs_dist, lqs_rq, lqs_time[1], lqs_time[2], lqs_time[3], lqs_lts, sep=",")
      write(lqs_out, file=paste(resloc, "/lqsi", i, ".csv", sep=""), ncol=length(lqs_out))
 
  } else {
      pca_time <- proc.time() 
      my_means <- colMeans(X)
      my_pca <- prcomp(X)
      pca_norm <- my_pca$rotation[,ncol(X)]
      pca_intercept <- -sum(my_means*pca_norm) # hyperplane equation is w^Tx + b = 0
      pca_dist <- sum(apply(X[1:m,], 1, hyper_dist_sq, w=pca_norm, b=pca_intercept))
      pca_res <- apply(X, 1, hyper_dist_sq, w=pca_norm, b=pca_intercept)
      pca_rq <- quantile(abs(pca_res), q)
      pca_lts <- sum(sort(pca_res)[1:m])
      pca_time <- proc.time() - pca_time

      pca_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "pca", pca_dist, pca_rq, pca_time[1], pca_time[2], pca_time[3], pca_lts, sep=",")
      print(pca_out)
      write(pca_out, file=paste(resloc, "/pcai", i, ".csv", sep=""))

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
          hbreg_time <- proc.time() - hbreg_time
          
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
          arob_times <- sapply(hbreg_results, function (x) x$arob_time)
          arob_time <- sum(arob_times)

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
          bb_times <- sapply(hbreg_results, function (x) x$bb_time)
          bb_time <- sum(bb_times)

          hbreg_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "hbreg", hbreg_dist, hbreg_rq, hbreg_time[1], hbreg_time[2], hbreg_time[3], hbreg_lts, sep=",")
          write(hbreg_out, file=paste(resloc, "/hbregi", i, ".csv", sep=""))
          arob_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "arob", arob_dist, arob_rq, arob_time[1], arob_time[2], arob_time[3], arob_lts, sep=",")
          write(arob_out, file=paste(resloc, "/arobi", i, ".csv", sep=""))
          bb_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "bb", bb_dist, bb_rq, bb_time[1], bb_time[2], bb_time[3], bb_lts, sep=",")
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

  # run CB and alg3
  if (dep_var == TRUE) {
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",",m, ",true,'", resloc, "')", sep="")) # dep_var = TRUE
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "','PCA')", sep=""))
  } else {
    print(add_path)
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",",m, ",false,'", resloc, "')", sep="")) # dep_var = TRUE
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','PCA')", sep=""))
  }


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
  

  q <- floor(q*nrow(X))
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="")
  print(paste("addpath('", srcloc, "');", sep=""))

  # CB and alg3
  if (dep_var == TRUE) {
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "')", sep=""))
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", resloc, "','PCA')", sep=""))
  } else {
    run_matlab_code(paste(add_path, " ", "run_cb(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "')", sep=""))
    run_matlab_code(paste(add_path, " ", "run_alg3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",false,'", resloc, "','PCA')", sep=""))
  }

}

## -----------------------------------------------------------------------------------------------------
get_mh <- function(dataloc, srcloc, fname, q, dep_var, timelimit, resloc) {
  #loc <- 'vary_everything'
  #fname <- 'm102n47m_outliers5num_clust1same_sideFALSEoutlier_dist1000i848.csv'
  #fname <- "m650n10m_outliers350i0.csv"
  #loc <- 'rd'
  #library(MASS)
  #source("../2021_02_hyper_data/olive/mpack.txt")
  #q <- 0.64
  print("mh")
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
      # M-estimator Huber
      mh_time <- proc.time()
      my_mh <- rlm(V1 ~ ., data=X, method="M")  
      mh_dist <- sum(my_mh$residuals[1:m]^2)
      mh_rq <- quantile(abs(my_mh$residuals), q)
      mh_lts <- sum(sort(my_mh$residuals^2)[1:m])
      mh_time <- proc.time() - mh_time

      # update after here
      mh_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "mh", mh_dist, mh_rq, mh_time[1], mh_time[2], mh_time[3], mh_lts, sep=",")
      print(mh_out)
      write(mh_out, file=paste(resloc, "/mhi", i, ".csv", sep=""), ncol=length(mh_out))
  } else {
      mh_time <- proc.time()
      mh_results <- lapply(1:n, fit_mh, X, m, n, q)
      mh_dists <- sapply(mh_results, function(x) x$mh_dist)
      mh_dist <- min(mh_dists)
      mh_norms <- sapply(mh_results, function(x) x$mh_norm)
      mh_intercepts <- sapply(mh_results, function(x) x$mh_intercept)
      mh_norm <- mh_norms[,which(mh_dists == min(mh_dists))]
      mh_intercept <- mh_intercepts[which(mh_dists == min(mh_dists))]
      mh_rqs <- sapply(mh_results, function(x) x$mh_rq)
      mh_rq <- mh_rqs[which(mh_dists == min(mh_dists))]
      mh_ltss <- sapply(mh_results, function(x) x$mh_lts)
      mh_lts <- mh_ltss[which(mh_dists == min(mh_dists))]
      mh_time <- proc.time() - mh_time

      mh_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "mh", mh_dist, mh_rq, mh_time[1], mh_time[2], mh_time[3], mh_lts, sep=",")
      print(mh_out)
      write(mh_out, file=paste(resloc, "/mhi", i, ".csv", sep=""))

  }
}

get_mm <- function(dataloc, srcloc, fname, q, dep_var, timelimit, resloc) {
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
      # MM-estimator
      mm_time <- proc.time()
      my_mm <- lmRob(V1 ~ ., data=X)  
      mm_dist <- sum(my_mm$residuals[1:m]^2)
      mm_rq <- quantile(abs(my_mm$residuals), q)
      mm_lts <- sum(sort(my_mm$residuals^2)[1:m])
      mm_time <- proc.time() - mm_time


      print(mm_out)
      write(mm_out, file=paste(resloc, "/mmi", i, ".csv", sep=""), ncol=length(mm_out))
  } else {
      mm_time <- proc.time()
      mm_results <- lapply(1:n, fit_mm, X, m, n, q)
      mm_dists <- sapply(mm_results, function(x) x$mm_dist)
      mm_dist <- min(mm_dists)
      mm_norms <- sapply(mm_results, function(x) x$mm_norm)
      mm_intercepts <- sapply(mm_results, function(x) x$mm_intercept)
      mm_norm <- mm_norms[,which(mm_dists == min(mm_dists))]
      mm_intercept <- mm_intercepts[which(mm_dists == min(mm_dists))]
      mm_rqs <- sapply(mm_results, function(x) x$mm_rq)
      mm_rq <- mm_rqs[which(mm_dists == min(mm_dists))]
      mm_ltss <- sapply(mm_results, function(x) x$mm_lts)
      mm_lts <- mm_ltss[which(mm_dists == min(mm_dists))]
      mm_time <- proc.time() - mm_time

      mm_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "mm", mm_dist, mm_rq, mm_time[1], mm_time[2], mm_time[3], mm_lts, sep=",")
      print(mm_out)
      write(mm_out, file=paste(resloc, "/mmi", i, ".csv", sep=""))
  }
}

get_rewlse <- function(dataloc, srcloc, fname, q, dep_var, timelimit, resloc) {
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
      # REWLSE-estimator
      rewlse_time <- proc.time()
      my_rewlse <- lmRob(V1 ~ ., data=X, control=lmRob.control(final.alg="Adaptive"))  
      rewlse_dist <- sum(my_rewlse$residuals[1:m]^2)
      rewlse_rq <- quantile(abs(my_rewlse$residuals), q)
      rewlse_lts <- sum(sort(my_rewlse$residuals^2)[1:m])
      rewlse_time <- proc.time() - rewlse_time

      rewlse_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "rewlse", rewlse_dist, rewlse_rq, rewlse_time[1], rewlse_time[2], rewlse_time[3], rewlse_lts, sep=",")
      print(rewlse_out)
      write(rewlse_out, file=paste(resloc, "/rewlsei", i, ".csv", sep=""), ncol=length(rewlse_out))
  } else {
      rewlse_time <- proc.time()
      rewlse_results <- lapply(1:n, fit_rewlse, X, m, n, q)
      rewlse_dists <- sapply(rewlse_results, function(x) x$rewlse_dist)
      rewlse_dist <- min(rewlse_dists)
      rewlse_norms <- sapply(rewlse_results, function(x) x$rewlse_norm)
      rewlse_intercepts <- sapply(rewlse_results, function(x) x$rewlse_intercept)
      rewlse_norm <- rewlse_norms[,which(rewlse_dists == min(rewlse_dists))]
      rewlse_intercept <- rewlse_intercepts[which(rewlse_dists == min(rewlse_dists))]
      rewlse_rqs <- sapply(rewlse_results, function(x) x$rewlse_rq)
      rewlse_rq <- rewlse_rqs[which(rewlse_dists == min(rewlse_dists))]
      rewlse_ltss <- sapply(rewlse_results, function(x) x$rewlse_lts)
      rewlse_lts <- rewlse_ltss[which(rewlse_dists == min(rewlse_dists))]
      rewlse_time <- proc.time() - rewlse_time

      rewlse_out <- paste(paste(dataloc, "/", fname, sep=""), i, nrow(X), n, m, q, "rewlse", rewlse_dist, rewlse_rq, rewlse_time[1], rewlse_time[2], rewlse_time[3], rewlse_lts, sep=",")
      print(rewlse_out)
      write(rewlse_out, file=paste(resloc, "/rewlsei", i, ".csv", sep=""))
  }
}

## -----------------------------------------------------------------------------------------------------
hyper_dist_sq <- function(w,b,x) {
  (sum(w*x) + b)^2/sum(w*w)
}



