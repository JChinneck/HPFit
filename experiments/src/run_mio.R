# The main function is run_mio() that calls one of three methods in
# MATLAB for minimizing the q^th absolute residual.  The methods are:
# MIO-BM from Bertsimas and Mazumder (2014), MIO1, an improved
# formulation developed by JWC, and MIO3, a three-phase approach 
# that uses MIO1, PCA, and re-instatement of points.  MIO-BM and 
# MIO1 are implemented in mio.m.  MIO3 is in mio3.m.  All three use
# mixed-integer optimization.  
#
# Main version and options:
# dep_var:
#   - TRUE: a dependent variable is specified, as in ordinary 
#           regression.  A residual is measured as the absolute 
#           difference between the reponse value and the 
#           vertical projection of the point on the fitted hyperplane
#   - FALSE: no dependent is specified. The intercept term is 
#            fixed to n, the original number of variables in the
#            dataset.  In the end, a residual is measured as the  
#            distance of a point to its orthogonal projection.
# q: 
#   - the percentile used to evaluate the fit of a hyperplane in
#     in the MIO1 and MIO-BM and phase 1 MIO3 formulations
#   - for CB-MIO3, q = 0 will use the outFinder q
#   - for CB-MIO3, q = "qout" will use the number of outliers on CB 
#     exit
# formulation:
#   - mio-bm: use the model proposed by Bertsimas and Mazumder (2014)
#   - mio1: use JWC's compact MIO model.
#   - mio3: use the three phase approach.
#
# calc_lqs_beta: 
#   - TRUE: use R's implementation of LQS to generate a warmstart
#           solution for MIO.  If dep_var=TRUE, allow each variable
#           to serve as the response in turn.
#   - FALSE: do not generate a warmstart using LQS.
#
# NOTES:
# - when a dep_var=TRUE, the response coefficient is -1.  When one is
#   not specified, the coefficients are normalized so that the 
#   intercept is n (beta_0 + beta^T x = 0).
# - mio.m and mio3.m create output files
#
# INPUTS:
# - options mentioned above
# - dataloc: folder containing data file
# - srcloc: folder containing mio.m and mio3.m
# - fname: data file name
# - timelimit: timelimit for the MIP solver
# - resloc: folder for the output 


# function for calculating the squared distance of a point to a
# hyperplane
hyper_dist_sq <- function(w,b,x) { 
  (sum(w*x) + b)^2/sum(w*w)
}


# function for applying LQS to a dataset with the j^th column 
# serving as the response.
fit_lqs <- function(j, X, m, n, q) {
  cat(j, " ")
  my_formula <- as.formula(paste("V", j, " ~ .", sep="")) # Vj is the response; all other variables are predictors
  if (q==1.0) {
    q=(nrow(X)-1)/nrow(X) 
  }
  my_lqs <- lqs(my_formula, data=X, method="lqs", quantile=floor(q*nrow(X)))  # Fit LQS model
  lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)] # Get coefficients excluding the intercept
  lqs_norm[is.na(lqs_norm)] <- 0.0
  lqs_norm <- append(lqs_norm, -1.0, after=j-1) # The coefficient for the response is -1
  lqs_intercept <- my_lqs$coefficients[1] # get the intercept
  lqs_dist <- sum(apply(X[1:m,],
                       1,
                       hyper_dist_sq,
                       w=lqs_norm,
                       b=lqs_intercept)) # sum of squared distances to the hyperplane for non-outliers
  lqs_res <- apply(X,1,hyper_dist_sq,w=lqs_norm,b=lqs_intercept) # sum of squared distances for all points
  lqs_rq <- quantile(lqs_res, q) # get the q^th percentile residual
  return(
         list(
              lqs_dist=lqs_dist, # non-outlier distances
              lqs_norm=lqs_norm, # coefficients
              lqs_intercept=lqs_intercept, # intercept
              lqs_rq=lqs_rq # q^th percentile residual
              )
  )
}

  
run_mio <- function(dataloc, srcloc, fname, q, dep_var, formulation, timelimit, resloc, calc_lqs_beta, mosekloc) {
  my_regexec <- regexec("m([0-9]+)n([0-9]+).+i([0-9]+).+csv", fname) # get m, n, i from the dataset name
  my_regmatch <- regmatches(fname, my_regexec)
  m <- as.numeric(my_regmatch[[1]][2])
  n <- as.numeric(my_regmatch[[1]][3])
  i <- as.numeric(my_regmatch[[1]][4])
  cat("\n", i, "\n")
  X <- read.csv(paste(dataloc,"/",fname,sep=""), header=FALSE) # read data

  lqs_beta <- -10
  if (calc_lqs_beta == TRUE) {
      if (dep_var == TRUE) { # response is first variable
          if (q==1.0) {
            q<-nrow(X)-1
          }
          my_lqs <- lqs(V1 ~ ., data=X, method="lqs", quantile=floor(q*nrow(X))) # fit LQS model for response; the first variable
          lqs_norm <- my_lqs$coefficients[2:length(my_lqs$coefficients)] # get coefficients but not intercept
          lqs_norm[is.na(lqs_norm)] <- 0.0
          lqs_intercept <- my_lqs$coefficients[1] # get intercept
          lqs_beta <- matrix(c(-1.0, lqs_intercept, lqs_norm),
                             nrow=length(lqs_norm)+2,
                             ncol=1) # beta_1 = -1 - the response
      } else {
          lqs_dist_min <- 10^8
          lqs_results <- lapply(1:n, fit_lqs, X, m, n, q)  # apply LQS with each variable 1...n as the response
          lqs_dists <- sapply(lqs_results, function(x) x$lqs_dist) # get all the distances
          lqs_dist <- min(lqs_dists) # get the best one

          lqs_norms <- sapply(lqs_results, function(x) x$lqs_norm) # get all the coefficients
          lqs_intercepts <- sapply(lqs_results, function(x) x$lqs_intercept) # get the intercepts
          lqs_norm <- lqs_norms[,which(lqs_dists == min(lqs_dists))] # get the coefficients of the best hyperplane
          lqs_intercept <- lqs_intercepts[which(lqs_dists == min(lqs_dists))] # get the intercept of the best hyperplane
          lqs_beta <- matrix(c(lqs_intercept, lqs_norm),
                             nrow=length(lqs_norm)+1,
                             ncol=1) # put intercept and coefficients together
          lqs_beta <- (lqs_beta/lqs_beta[1,1])*n  # normalize so that beta_0 = n
      }
  }

  if (class(q) == "numeric") {
    q <- floor(q*nrow(X)) # Convert q from percentile to order statistic for MATLAB implementations
  }
  make_lqs_beta <- rmat_to_matlab_mat(lqs_beta, matname="lqs_beta") # create the warm start from LQS
  #add_path <- paste("addpath('", srcloc, "');", sep="") 
  add_path <- paste("addpath('", srcloc, "','", mosekloc, "');", sep="") # add the path to the MATLAB files and MOSEK
  if (formulation == "mio3" | formulation == "lqs-mio3" | formulation == "alg3-mio3") { # three-phase approach
    if (dep_var == TRUE) { # first variable is response
      cat("running mio3 ")
      run_matlab_code(paste(add_path, " ", make_lqs_beta, " ", "mio3(", i, ",'",dataloc, "/", fname,"',",q,",lqs_beta,", m, ",true,'", formulation, "','", resloc, "',", timelimit, ",-1)", sep="")) # dep_var = TRUE
    } else { 
      run_matlab_code(paste(add_path, " ", make_lqs_beta, " ", "mio3(", i, ",'",dataloc, "/", fname,"',",q,",lqs_beta,", m, ",false,'", formulation, "','", resloc, "',", timelimit, ",-1)", sep="")) # dep_var = FALSE
    }
     
  } else if (formulation == "cbmio3") {
    if (dep_var == TRUE) { # first variable is response
      cat("running cbmio3 ")
      cat(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"',", q,",", m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
      print(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"',",q,",", m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
      if (class(q) == "numeric") {
        run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"',",-q,",",m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
      } else {
        run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"','qout',",m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
      }
    } else {
      if (class(q) == "numeric") {
        run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"',",-q,",", m, ",false,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = FALSE
      } else {
        run_matlab_code(paste(add_path, " ", "run_cbmio3(", i, ",'",dataloc, "/", fname,"','qout',", m, ",false,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = FALSE
      }
    }
  } else{ # MIO1 or MIO-BM
    if (dep_var == TRUE) { # first variable is response 
      print(paste(add_path, " ", make_lqs_beta, " ", "mio(", i, ",'",dataloc, "/", fname,"',",q,",lqs_beta,", m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) 
      run_matlab_code(paste(add_path, " ", make_lqs_beta, " ", "mio(", i, ",'",dataloc, "/", fname,"',",q,",lqs_beta,", m, ",true,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = TRUE
    } else {
      run_matlab_code(paste(add_path, " ", make_lqs_beta, " ", "mio(", i, ",'",dataloc, "/", fname,"',",q,",lqs_beta,", m, ",false,'", formulation, "','", resloc, "',", timelimit, ")", sep="")) # dep_var = FALSE
    }
  }
}


