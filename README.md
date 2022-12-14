# HPFit
Matlab code for fitting hyperplanes (for regression, as well as for general data sets). The heuristic methods are fast and effective for very large data sets. The mixed-integer optimizations are slow but exact.

A vital step is the identification and removal of outlier points. All of the following codes rely on a robust method for identifying outliers.

**CBgen** fits a hyperplane to a general set of data, formatted in a table with rows for points and columns for variables. The goal is to find a fit that is as close to as many points as possible.

**CBreg** fits a regression hyperplane so that it is as close to as many of the output variable values as possible.

**CBqgen** fits a hyperplane that minimizes the qth largest error in a data set, where error is the Euclidean distance from the hyperplane. q can be expressed as a percentile (often the median).

**CBqreg** fits a regression hyperplane that minimizes the qth largest residual. q can be expressed as a percentile (often the median).

**CBMIOgen** is a mixed-integer optimization that finds a hyperplane to exactly minimize the qth largest error in a data set, where error is the Euclidean distance from the hyperplane. q can be expressed as a percentile (often the median).

**CBMIOreg** is a mixed-integer optimization that finds a regression hyperplane to exactly minimize the qth largest residual. q can be expressed as a percentile (often the median).
