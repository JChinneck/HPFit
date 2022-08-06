% run_cbmio3() calls one of CBMIOreg.m or CBMIOgen.m depending on
% whether there is a dependent variable or not.
%
% Main versions and options:
% dep_var:
%   - true: a dependent variable is specified, as in ordinary 
%           regression.  A residual is measured as the squared 
%           difference between the reponse value and the 
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent is specified. The intercept term is 
%            fixed to n, the original number of variables in the
%            dataset.  In the end, a residual is measured as the  
%            distance of a point to its orthogonal projection.
% q: 
%   - the order statistic used to evaluate the fit of a hyperplane
%     in Phase 1.  The fit is determined by the q^th largest residual.
% NOTES:
% - when a dependent variable is specified, the corresponding 
%   coefficient is set to -1.  When one is not specified, the 
%   intercept is set to n.
% - a result datafile is created containing the data file with 
% bnd2,  gurobi runtime, gamma, 
% TSEstar after CB, TSE after CB, gamma after CB
% TSEstar after MIO, TSE after MIO, gamma after MIO
% TSEstar after PCA on MIO, TSE after PCA on MIO, gamma after PCA on MIO
% TSEstar after PCA on MIO, TSE after PCA on MIO, gamma after PCA on MIO
% TSEstar after final PCA, TSE after final PCA, gamma after final PCA
% MIO time,  CB time, total time
%
% INPUTS:
% - options mentioned above
% - iteration: iteration number.  Used in the output filename.
% - datafname: full path to data file.
%   for MIO3
% - m_normal: the  number of non-outlier rows of data.
% - resloc: path to folder where output file will reside.
% - timelimit: time limit for MIP solver.
%

function  run_cbqmio(iteration, datafname, q, m_normal, dep_var, formulation, resloc, timelimit)
disp("running cbqmio")
dep_var

X = readtable(datafname);
X = X{:,:}; % convert data table to an array
[m,n] = size(X);  % the original size of the dataset

gbparams.TimeLimit = timelimit;
qparams.q = -q;

out_fname = strcat(resloc, "/", formulation, "i", int2str(iteration), ".csv")
disp(out_fname)
out_file = fopen(out_fname, 'w');

if dep_var == true
    qparams.maxResid = -16;
    y = X(:,1); % first column is response
    A = X(:,2:n);
    [output, inc] = CBqreg(y, A, qparams)
    output.w
    cbq_beta = [-1.0;  output.w0; output.w];
    cbq_beta
else 
    qparams.maxDist = -16;
    [output, inc] = CBqgen(X, qparams)
    cbq_beta = [-output.RHS; output.weights];
    cbq_beta = (cbq_beta/cbq_beta(1,1))*n
end
[beta_star, f_beta_star] =mio(iteration, datafname,q, cbq_beta, m_normal, dep_var, formulation, resloc, timelimit);


fclose(out_file);

return
end
