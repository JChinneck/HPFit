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

function  run_cbmio3(iteration, datafname, q, m_normal, dep_var, formulation, resloc, timelimit)
disp("running cbmio3")
dep_var

X = readtable(datafname);
X = X{:,:}; % convert data table to an array
[m,n] = size(X);  % the original size of the dataset

gbparams.TimeLimit = timelimit;
mioparams.mtru = m_normal;
mioparams.q = q;
mioparams.maxResid = -16;
mioparams.maxDist = -16;

out_fname = strcat(resloc, "/", formulation, "i", int2str(iteration), ".csv")
disp(out_fname)
out_file = fopen(out_fname, 'w');

if dep_var == true
    y = X(:,1); % first column is response
    A = X(:,2:n);
    [result, output] = CBMIOreg(y, A, mioparams, gbparams)
    q = output.q; 
% filename, iteration, total number of points, number of variables, number of non-outliers, percentile for LQS, formulation- cbmio3, 
% gurobi runtime, gurobi status, CB gamma, objbound
% TSEstar after CB, TSE after CB, gamma after CB
% TSEstar after MIO, TSE after MIO, gamma after MIO
% TSEstar after PCA on MIO, TSE after PCA on MIO, gamma after PCA on MIO
% TSEstar after PCA on MIO, TSE after PCA on MIO, gamma after PCA on MIO
% TSEstar after final PCA, TSE after final PCA, gamma after final PCA
% MIO time,  CB time, total time
    fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%f,%f,%f,%d,%d\n", datafname, iteration, m, n, m_normal, q, formulation, result.runtime, result.status, output.gamma(4,1), result.objbound, output.TSEstar(1,1), output.TSE(1,1), output.gamma(1,1), output.TSEstar(2,1), output.TSE(2,1), output.gamma(2,1), output.TSEstar(3,1), output.TSE(3,1), output.gamma(3,1), output.TSEstar(4,1), output.TSE(4,1), output.gamma(4,1),output.closeAll(1,1),output.closeAll(2,1),output.closeAll(3,1),output.closeAll(4,1),output.MIOTime, output.CBregTime, output.solTime,output.cbq,output.outfinderq); 
else 
    [result, output] = CBMIOgen(X, mioparams, gbparams)
    q = output.q;
    fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%f,%f,%f,%d,%d\n", datafname, iteration, m, n, m_normal, q, formulation, result.runtime, result.status, output.gamma(4,1), result.objbound, output.TSEstar(1,1), output.TSE(1,1), output.gamma(1,1), output.TSEstar(2,1), output.TSE(2,1), output.gamma(2,1), output.TSEstar(3,1), output.TSE(3,1), output.gamma(3,1), output.TSEstar(4,1), output.TSE(4,1), output.gamma(4,1),output.closeAll(1,1),output.closeAll(2,1),output.closeAll(3,1),output.closeAll(4,1),output.MIOTime, output.CBgenTime, output.solTime, output.cbq,output.outfinderq ); 
end


fclose(out_file);

return
end
