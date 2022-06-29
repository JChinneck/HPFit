% mio() uses MIP to fit a hyperplane that minimizes the q^th residual.
% Purpose: to find the best fitting hyperplance to a set of input 
% data, where best fitting means that the q^th residual is minimized.
% which can be affected by outliers.
% Two general formulations are presented.  One is based on the 
% original presented in Bertsimas and Mazumder 2014 (2.11). The other
% is an equivalent but more compact formulation developed by JC. For 
% each formulation, there are two versions: 1) a dependent variable is
% specified (like the original published version) and 2) no dependent
% variable is specified.  
%
% Main versions and options:
% formulation:
%   - mio-bm: the original MIP formulation presented in (2.11) of 
%             Bertsimas and Mazumder (2014).  Warm start with alg3.
%   - mio1: an equivalent but more compact formualtion developed by
%           JC.  Warm start with alg3.
%   - lqs-mio-bm: mio-bm warmstarted with LQS.
%   - lqs-mio1: mio1 warmstarted with LQS.
% dep_var:
%   - true: a dependent variable is specified, as in ordinary 
%           regression.  A residual is measured as the absolute 
%           difference between the reponse value and the 
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent is specified. The intercept term is 
%            fixed to n, the original number of variables in the
%            dataset.  A residual is measured as the  
%            distance of a point to its orthogonal projection.
% q: 
%   - the order statistic used to evaluate the fit of a hyperplane. 
%     The fit is determined by the q^th largest residual.
% NOTES:
% - when a dependent variable is specified, the corresponding 
%   coefficient is set to -1.  When one is not specified, the 
%   intercept is set to n.
% - a result datafile is created containing the data file with full 
%   path, iteration number, number of rows, number of variables, 
%   number of non-outliers, q, formulation, total squared error, 
%   MIP runtime, MIP status, gamma, best MIP bound, number of outliers
%   identified as one of the q smallest by best MIP feasible solution.
% 
% INPUTS:
% - options mentioned above
% - iteration: iteration number.  Used in the output filename.
% - datafname: full path to data file.
% - lqs_beta: an initial solution generated using LQS as implemented
%             in R. If -10, not lqs_beta provided.
% - m_normal: number of non-outlier rows of data.  After that they are
%             outliers.
% - resloc: path to folder where output file will reside.
% - timelimit: time limit for MIP solver.
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
 

function [beta_star, f_beta_star] = mio(iteration, datafname,q, lqs_beta, m_normal, dep_var, formulation, resloc, timelimit)

X = readtable(datafname);
X = X{:,:};
[m,n] = size(X); 
if dep_var == true
    X = [X(:,1) ones(m,1) X(:,2:n)]; % add column of 1s; first column is response, second column is 1s to represent intercept, then 2:n are data
else
    X = [ones(m,1) X]; % add column of 1s at beginning for the intercept, then 2:n are data
end
[m,n] = size(X); % now n will include the intercept column and is one more than the original data


% Set up the MIP
if strcmp(formulation, "mio-bm") | strcmp(formulation, "lqs-mio-bm")
    model.obj = [1.0; zeros(2*m,1) ; zeros(m,1); zeros(m,1); zeros(m,1); zeros(n,1)]; % gamma ; rplus/rminus; mu; mubar; z; beta
    model.lb  = [zeros(1+5*m,1); -inf(n,1)];
    if dep_var == true % dependent variable - regression; first variable is response
        model.A   = [ 
                  sparse(-ones(m,1)) speye(m) speye(m) -speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n))  ; % rplus + rminus - gamma = mubar - mu, for each point i
                  sparse(zeros(m,1)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(X) ; % rplus - rminus = y-x^T beta, for each point i
                  sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)) ; % sum of zs is q
                  sparse(ones(m,1)) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(zeros(m,m)) -speye(m) sparse(zeros(m,m)) sparse(zeros(m,n))  ; % mu <= gamma, for each point i
                  sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1))] ; % beta_1=-1
      
        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; zeros(m,1) ; -1];
        %model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; zeros(m,1) ; ones(m,1) ; -1];
    else % no dependent variable
        model.A   = [ 
                  sparse(-ones(m,1)) speye(m) speye(m) -speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % rplus + rminus - gamma = mubar - mu, for each point i
                  sparse(zeros(m,1)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(X) ; % rplus - rminus = x^T beta, for each point i
                  sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)) ; % sum of zs is q
                  sparse(ones(m,1)) sparse(zeros(m,m)) sparse(zeros(m,m)) sparse(zeros(m,m)) -speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % mu <= gamma for each point i
                  sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)) ] ; % beta_0=n; n-1 because we have a column of 1s

        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; zeros(m,1) ; n-1];
    end
    model.sense = [repmat('=', m, 1); repmat('=', m, 1); '='; repmat('>', m, 1) ; '='] ;
    model.vtype = ['C'; repmat('C', 4*m, 1) ; repmat('B', m, 1); repmat('C', n, 1) ];
    model.varnames = cellstr(['gamma' ; repmat('rplus',m, 1) + string(1:m)' ; repmat('rminus',m, 1) + string(1:m)' ; repmat('mu',m, 1) + string(1:m)' ; repmat('mubar',m, 1) + string(1:m)' ; repmat('z',m, 1) + string(1:m)' ; repmat('beta',n, 1) + string(1:n)']) 
else % formulation is MIO1
    model.obj = [1.0; zeros(m,1) ; zeros(m,1); zeros(m,1); zeros(m,1); zeros(n,1)]; % gamma ; r; eplus; eminus; z; beta
    model.lb  = [zeros(1+4*m,1); -inf(n,1)];
    if dep_var == true % dependent variable - regression; first variable is response
        model.A   = [ sparse(zeros(m,1)) sparse(zeros(m,m)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(X) ; % beta^T x - y + eplus - eminus = 0, for each point i
                      sparse(-ones(m,1)) -speye(m) speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % eplus + eminus - r - gamma <= 0, for each point i
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)); % sum of zs is q
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)) ] ; % beta_1=-1
        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; -1 ];
    else
        model.A   = [ sparse(zeros(m,1)) sparse(zeros(m,m)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(X) ; % beta^Tx + eplus - eminus = 0, for each point i
                      sparse(-ones(m,1)) -speye(m) speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % eplus + eminus - r - gamma <= 0, for each point i
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)); % sum of zs is q
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)) ] ; % beta_0=n
        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; n-1 ]; % n-1 because we have a column of 1s for the intercept
    end
    model.sense = [repmat('=', m, 1); repmat('<', m, 1); '='; '='] ;
    model.vtype = ['C'; repmat('C', 3*m, 1) ; repmat('B', m, 1); repmat('C', n, 1)];
    model.varnames = cellstr(['gamma' ; repmat('r',m, 1) + string(1:m)' ; repmat('eplus',m, 1) + string(1:m)' ; repmat('eminus',m, 1) + string(1:m)' ; repmat('z',m, 1) + string(1:m)' ; repmat('beta',n, 1) + string(1:n)']) 
end
for k=1:m
    if strcmp(formulation, "mio-bm") | strcmp(formulation, "lqs-mio-bm")
        model.sos(k).type = 1;
        model.sos(k).index = [(1+2*m+m+k) (1+2*m+k)]'; % mubar, mu
        model.sos(m+k).type = 1;
        model.sos(m+k).index = [(1+k) (1+m+k)]'; % rplus, rminus
        model.sos(2*m+k).type = 1;
        model.sos(2*m+k).index = [(1+2*m+2*m+k) (1+2*m+k)]'; % z, mu
    else % MIO1
        model.sos(k).type = 1;
        model.sos(k).index = [(1+k) (1+3*m+k)]'; % r, z
    end
end
disp("alg 3 start")
[beta1, f_beta1] = algorithm3(X, q, dep_var, "PCA"); % algorithm 3 is given by Bertsimas and Mazumder as a way to derive initial solutions
beta1
disp("alg 3 end")
if strcmp(formulation, "mio-bm")
    model.StartNumber = 0;
    model.start = [nan; NaN(5*m,1);   beta1];  % from algorithm 3
else %MIO1
    model.StartNumber = 0;
    model.start = [nan; NaN(4*m,1);   beta1];  % from algorithm 3

end

if strcmp(formulation, "lqs-mio-bm")% if a solution given by R's LQS method is given
    model.StartNumber = 1;
    model.start = [nan; NaN(5*m,1) ; lqs_beta]; 
end
if strcmp(formulation, "lqs-mio1")
    model.StartNumber = 1;
    model.start = [nan; NaN(4*m,1) ; lqs_beta]; 
end

model.modelsense = 'min';
%gurobi_write(model, 'mio.lp');
params = struct();
%params.OutputFlag = 0;
if timelimit > 0
    params.TimeLimit = timelimit;
end
params.Symmetry = 2;
params.Threads = 1
disp("solving")
result = gurobi(model, params);
result.status
if strcmp(result.status, 'OPTIMAL')
    if strcmp(formulation, "mio-bm") | strcmp(formulation, "lqs-mio-bm")
        beta_star = result.x((1+5*m+1):(1+5*m+n),1);
        z = result.x(1+4*m+1:1+4*m+m,1)
    else % MIO1
        beta_star = result.x((1+4*m+1):(1+4*m+n),1);
        z = result.x(1+3*m+1:1+3*m+m,1)
    end
else 
    if result.mipgap ~= Inf
        fprintf("Using incumbent solution\n")
        if strcmp(formulation, "mio-bm") | strcmp(formulation, "lqs-mio-bm")
            beta_star = result.pool(1).xn((1+5*m+1):(1+5*m+n),1);
            z = result.pool(1).xn(1+4*m+1:1+4*m+m,1)
        else % MIO1
            beta_star = result.pool(1).xn((1+4*m+1):(1+4*m+n),1);
            z = result.pool(1).xn(1+3*m+1:1+3*m+m,1)
        end
    else
        fprintf("No solution available\n")
        return
    end
end
f_beta_star = result.objval;
qth_residuals = [f_beta1 f_beta_star];

num_outliers_in_q = sum(z((m_normal+1):m,1))

% get sum of squared error on non outliers
dist = abs(X*beta_star); 
if dep_var == true % get error along response direction, recall that beta_1 = -1
    tot_err = sum(dist(1:m_normal,1).*dist(1:m_normal,1));
    sorteddist = sort(dist(:,1).*dist(:,1));
    tsestar = sum(sorteddist(1:m_normal))
    tse = sum(sorteddist(1:q))
else % get orthogonal error
    gradLen = norm(beta_star(2:n,1)); % first coefficient is the intercept; exclude that from the gradLen calculation
    edist = abs(dist/gradLen);
    sortededist = sort(edist(:,1).*edist(:,1));
    tsestar = sum(sortededist(1:m_normal))
    tse = sum(sortededist(1:q))
    tot_err = sum(edist(1:m_normal,1).*edist(1:m_normal,1));
end

out_fname = strcat(resloc, "/", formulation,"i",int2str(iteration), ".csv")
disp(out_fname)
out_file = fopen(out_fname, "w");

%disp(result.status)
beta_star
 
% filename, total number of points, number of variables, number of non-outliers, q - percentile for LQS, formulation - mio-bm or mio1,total squared error to hyperplane (along response or orthogonal, gurobi runtime, gamma 
fprintf(out_file, "%s,%d,%d,%d,%d,%d,%s,%f,%f,%s,%f,%f,%d,%f,%f\n", datafname, iteration, m, n-1, m_normal, q, formulation, tot_err, result.runtime, result.status, f_beta_star, result.objbound, num_outliers_in_q,tse,tsestar);

fclose(out_file);
return
        
end
