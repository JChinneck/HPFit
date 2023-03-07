% run_cbq() calls CBqreg.m or CBqgen.m for fitting hyperplanes
% Main options: 
% q: 
%   - the order statistic used to evaluate the fit of a hyperplane. 
%     The fit is determined by the q^th largest residual.
% dep_var:
%   - true: a dependent variable is specified, as in ordinary
%           regression.  A residual is measured as the squared
%           difference between the reponse value and the 
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent variable is specified. 
%            A residual is measured using the elastic LP criterion,
%            which is the absolute difference from zero of the inner
%            product of coefficients and data.  Performance at the 
%            end is measured using distance to the orthogonal 
%            projection.
% 
% NOTES:
% - when a dependent variable is specified, the corresponding 
%   coefficient is set to -1.  When one is not specified, the 
%   intercept is set to n.
%
% INPUTS: 
% - options mentioned above
% - iteration: iteration number for the dataset.  Used for the output
%              filename
% - datafname: full path to the data file name.
% - m_normal: the number of non-outliers in the data.
% - resloc: directory for the output file.
% - 
% 
% OUTPUTS:
% - beta_star: the coefficients of the best-fit hyperplane.  When 
%              dep_var=true, the first coefficient is -1 for the 
%              dependent variable, the second coefficient is the 
%              intercept, and the remaining are the coefficients for 
%              other variables.  When dep_var=false, the first 
%              coefficient is n for the intercept and the remaining
%              are the coefficients for the other variables.
% - f_beta_star: gamma, the objective function value.  The absolute
%                value of the q^th largest residual.
% - results are written to a file.

function [beta_star, f_beta_star] = run_cbq(iteration, datafname, q, m_normal, dep_var, resloc) 

inParam.mgood = m_normal;
inParam.maxDist = -16;
inParam.maxResid = -16;
inParam.q = -q;


disp("running cbq")
dep_var

X = readtable(datafname);
X = X{:,:}; % convert data table to an array
[m,n] = size(X);  % get the size of dataset
out_fname = strcat(resloc, "/cbqi", int2str(iteration), ".csv")
out_file = fopen(out_fname, "w");
if dep_var == true
    [output,inc] = CBqreg(X(:,1), X(:,2:n), inParam)
    % filename including path, total number of points, number of variables, number of non-outliers, q - percentile for LQS, formulation - cb, gamma, total squared error to hyperplane (along response or orthogonal), runtime, LTS 
    fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%f,%f,%f,%d,%d\n", datafname, iteration, m, n, m_normal, q, "cbq", output.gamma, -1.0, output.CBregTime, -1.0, -1,-1);
else
    [output,inc] = CBqgen(X, inParam)
    % filename including path, total number of points, number of variables, number of non-outliers, q - percentile for LQS, formulation - cb, gamma, total squared error to hyperplane (along response or orthogonal), runtime, LTS 
    fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%f,%f,%f,%d,%d\n", datafname, iteration, m, n, m_normal, q, "cbq", output.gammaN, -1.0, output.CBgenTime, -1.0, -1, -1);
end
 
fclose(out_file);
return

        
end



