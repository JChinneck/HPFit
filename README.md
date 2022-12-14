# HPFit
Matlab code for fitting hyperplanes (for regression, as well as for general data sets).

A vital step is the identification and removal of outlier points.

**CBgen** fits a hyperplane to a general set of data, formatted in a table with rows for points and columns for variables. The goal is to find a fit that is as close to as many points as possible.

**CBreg** fits a regression hyperplane so that it is as close to as many of the output variable values as possible.

**CBqgen** fits a hyperplane that minimizes the qth largest error, where error is the Euclidean distance from the hyperplane. q can be expressed as a percentile (often the median).

**CBqreg** fits a regression hyperplane that minimizes the qth largest residual. q can be expressed as a percentile (often the median).

