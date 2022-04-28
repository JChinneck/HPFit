% rbm() is a calling function for the RBM method for fitting a 
% hyperplane.  This function is called by run_competitors.R. RBM was
% developed by JWC; the "relative better measure" heuristic.  This 
% function cals RBMpcaV<#>() in RBMpcaV<#>.m where # is the version 
% number.  More details about RBM are there.
% q is the number of points for LQS
% m_normal is the number of non-outliers (at the top of the dataset)
% dep_var=true means there is a dependent variable 
%
% Main options:
% q: 
%   - the order statistic used to evaluate the fit of a hyperplane. 
%     The fit is determined by the q^th largest residual.  For 
%     comparisons to LTS/LQS methods.
% dep_var:
%   -true: a dependent variable is specified, as in ordinary 
%          regression.  A residual is measured as the absolute
%          difference between the response value and the
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent variable is specified. PCA is used as a 
%            subroutine to fit hyperplanes.  If insufficient 
%            data are available, an elastic LP is solved.
%            A residual is measured using the distance to the 
%            orthogonal projection.
% INPUTS: 
%   parameters mentioned above
% - iteration: iteration number for the dataset.  Used for the output
%              filename
% - datafname: full path to the data file name.
% - m_normal: the number of non-outliers in the data.
% - resloc: directory for the output file.
%
% OUTPUTS:
% - results are written to a file.  

function rbm(iteration, datafname, q, m_normal, dep_var, resloc)
formulation="rbm";

X = readtable(datafname);
X = X{:,:};
[m,n] = size(X); 
if dep_var == true
    X = [X(:,1) ones(m,1) X(:,2:n)]; % add column of 1s for the intercept; first column is response
else
    X = [ones(m,1) X]; % add column of 1s at beginning for intercept
end
[m,n] = size(X);  % n includes the intercept

RBMparam = struct(); % default settings from RBMpcaV1.m
RBMparam.feaTol=1.0e-6;
RBMparam.minFrac=0.5;
RBMparam.maxRemoveFrac=0.75;
RBMparam.maxDist=-16;
RBMparam.stopCondn='noBetter';
RBMparam.reinstate=1;


fprintf("running RBMpcaV1\n")
[errs,solTime,nHPs] = RBMpcaV1(X, m_normal, q, RBMparam,dep_var) % run RBM


% get sum of squared error on non outliers
if dep_var == true % get error along response direction
    beta_star = [-1.0; errs.weights(:,3)]; % add the response coefficient 
    dist = abs(X*beta_star); % response is first coefficient, intercept is second, regression coefficients follow
    tot_err = sum(dist(1:m_normal,1).*dist(1:m_normal,1)); % sum of squared distances along response direction
    dist_sort = sort(dist(:,1).*dist(:,1)); % sorted absolute distances
    rbm_lts = sum(dist_sort(1:m_normal)) % sum of m smallest absolute distances 
else
    tot_err = errs.totSqDistTru(1,3); % sum of squared orthogonal distances for true points
    dist = errs.edist(:,3).*errs.edist(:,3); % squared orthogonal distances
    dist_sort = sort(dist); % sorted squared distnces
    rbm_lts = sum(dist_sort(1:m_normal)) % sum of m smallest squared distances
    
end

out_fname = strcat(resloc, "/", formulation,"i",int2str(iteration), ".csv") % store results in rbmi<iteration>.csv
disp(out_fname)
out_file = fopen(out_fname, "w");


 
% filename, total number of points, number of variables, number of non-outliers, q - percentile for LQS, formulation - rbm,total squared error to hyperplane (along response or orthogonal), runtime, gamma 
fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%f,%f,%f\n", datafname, iteration, m, n-1, m_normal, q, formulation, tot_err, solTime, errs.gamma(1,3),rbm_lts);

fclose(out_file);
return

        
end




