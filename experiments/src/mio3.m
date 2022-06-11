% mio() uses a three-phase approach to fitting a hyperplane.  The 
% phases are:
% 1. Runs the improved MIO minimization of the qth percentile error, in the
% sorted set of errors, subject to the time limit.
% 2. Places a PCA hyperplane on the q points kept by the MIO solution.
% 3. Checks which points are outliers relative to the step 2 HP. Outliers
%    removed, nonoutliers kept, and a final HP is placed. Generally there
%    is likely to be only reintroduction of points that are not outliers.
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
% - a result datafile is created containing the data file with full 
%   path, iteration number, number of rows, number of variables, 
%   number of non-outliers, q, formulation, total squared error for 
%   non outliers, 
%   MIP runtime, MIP status, gamma, best MIP bound, total squared 
%   error for non-outliers after Phase 1, total squared error for all
%   points after Phase 1, total squared error for non-outliers after 
%   Phase 2, total squared error for all points after Phase 2, total
%   squared error for non-outlier points after Phase 3 (again), total
%   suared errrof for all points after Phase 3, number of outliers
%   identified as one of the q smallest by best MIP feasible solution,
%   algorithm 3 time.
%
% INPUTS:
% - options mentioned above
% - iteration: iteration number.  Used in the output filename.
% - datafname: full path to data file.
% - lqs_beta: initial solution generate using LQS in R; a warmstart 
%   for MIO3; if -10, none is provided
% - m_normal: the  number of non-outlier rows of data.
% - resloc: path to folder where output file will reside.
% - timelimit: time limit for MIP solver.
% - errs: object to store errors after each phase. On input, contains
%         maxDist; to count number of points within maxDist distance
% - formulation: "mio3" 
%
% OUTPUTS:
% - beta_star: the coefficients of the best-fit hyperplane after
%              phase 3.  When 
%              dep_var=true, the first coefficient is -1 for the 
%              dependent variable, the second coefficient is the 
%              intercept, and the remaining are the coefficients for 
%              other variables.  When dep_var=false, the first 
%              coefficient is n for the intercept and the remaining
%              are the coefficients for the other variables.
% - f_beta_star: gamma, the objective function value from the first 
%                phase.  The absolute value of the q^th largest residual.

function [beta_star, f_beta_star] = mio3(iteration, datafname,q, lqs_beta, m_normal, dep_var, formulation, resloc, timelimit, errs)
    disp("running mio3")

    X = readtable(datafname);
    X = X{:,:}; % convert data table to an array
    [m,n] = size(X);  % the original size of the dataset
    if dep_var == true
        X = [X(:,1) ones(m,1) X(:,2:n)]; % add column of 1s for the intercept after the response column
    else
        X = [ones(m,1) X]; % add column of 1s at beginning for the intercept
    end
    [m,n] = size(X);  % get the size with the intercept

    % STEP 1: set up and solve the MIO
    model.obj = [1.0; zeros(m,1) ; zeros(m,1); zeros(m,1); zeros(m,1); zeros(n,1)]; % gamma ; r; eplus; eminus; z; beta
    model.lb  = [zeros(1+4*m,1); -inf(n,1)];
    if dep_var == true % dependent variable - regression; first variable is response
        model.A   = [ sparse(zeros(m,1)) sparse(zeros(m,m)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(X) ; % beta^T x - y + eplus - eminus = 0, for each point i
                      sparse(-ones(m,1)) -speye(m) speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % eplus + eminus - r - gamma <= 0, for each point i
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)); % sum of zs is q
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)) ] ; % beta_1=-1
        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; -1 ];
    else
        model.A   = [ sparse(zeros(m,1)) sparse(zeros(m,m)) speye(m) -speye(m) sparse(zeros(m,m)) sparse(X) ; % Elastic cons, of form beta^T x + eplus - eminus = 0, for each point i 
                      sparse(-ones(m,1)) -speye(m) speye(m) speye(m) sparse(zeros(m,m)) sparse(zeros(m,n)) ; % eplus + eminus - r - gamma <= 0, for each point i
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(ones(1,m)) sparse(zeros(1,n)); % sum of zs is q
                      sparse(zeros(1,1)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)) ] ; % beta_0=n
        model.rhs = [ zeros(m,1) ; zeros(m,1) ; q ; n-1 ]; % n-1 because we have a column of 1s for the intercept
    end
    model.sense = [repmat('=', m, 1); repmat('<', m, 1); '='; '='] ;
    model.vtype = ['C'; repmat('C', 3*m, 1) ; repmat('B', m, 1); repmat('C', n, 1)];
    model.varnames = cellstr(['gamma' ; repmat('r',m, 1) + string(1:m)' ; repmat('eplus',m, 1) + string(1:m)' ; repmat('eminus',m, 1) + string(1:m)' ; repmat('z',m, 1) + string(1:m)' ; repmat('beta',n, 1) + string(1:n)']) 
    
    for k=1:m
         % mio1
        model.sos(k).type = 1;
        model.sos(k).index = [(1+k) (1+3*m+k)]'; % r,z
    end
 
    % get algorithm 3 warm start
    tic;
    if lqs_beta ~= -10 % only if an LQS warm start is provided
        disp("algorithm 3")
        [beta1, f_beta1] = algorithm3(X, q, dep_var, "PCA")
        disp("finished algorithm 3")
    end
    alg3_time = toc;
    %mio1
    lqs_beta
    if lqs_beta ~= -10
        model.StartNumber = 1; % from alg3
        model.start = [nan; NaN(4*m,1);   beta1];  % set alg3 warm start
        %model.StartNumber = 2; % from LQS
        %model.start = [nan; NaN(4*m,1) ; lqs_beta]; % set LQS warm start
    end

    model.modelsense = 'min';
    %gurobi_write(model, 'mio.lp');
    %disp("finished writing mio.lp")
    params = struct();
    %params.OutputFlag = 0;
    if timelimit > 0
        params.TimeLimit = timelimit;
    end
    params.Symmetry = 2;
    params.Threads = 1;
    disp("solving MIO")
    result = gurobi(model, params);
    disp("finished solving MIO")
    % mio1
    %beta_star = result.x((1+4*m+1):(1+4*m+n),1);
  

% About solutions:
% Check results.status:
%   If 'OPTIMAL' then take result.x
%   If 'TIME_LIMIT' then check whether result.mipgap is Inf
%     If yes, exit with no solution.
%     Else take result.pool.objval(1,1) as incumbent value of objective and
%       take result.pool(1).xn as the associated solution.
% get results
fprintf("\nSolution status: %s\n",result.status)
if strcmp(result.status, 'OPTIMAL')
    f_beta_star = result.objval;
    beta_star = result.x((1+4*m+1):(1+4*m+n),1);
%     eplus = result.x(n+1:n+m,1);
%     eminus = result.x(n+m+1:n+2*m,1);
%     rel = result.x(n+2*m+1:n+3*m,1);
    z = result.x((1+3*m+1):(1+3*m+m),1);
    output.gamma = result.x(1,1);
else
    if result.mipgap ~= Inf
        % Gurobi stopped for some reason, maybe time limit, but it has an
        % incumbent solution, so use that
        fprintf("Using incumbent solution\n")
        f_beta_star = result.objval;
        beta_star = result.pool(1).xn((1+4*m+1):(1+4*m+n),1);
%         eplus = result.pool(1).xn(n+1:n+m,1);
%         eminus = result.pool(1).xn(n+m+1:n+2*m,1);
%         rel = result.pool(1).xn(n+2*m+1:n+3*m,1);
        z = result.pool(1).xn(1+3*m+1:1+3*m+m,1);
        output.gamma = result.pool(1).xn(1,1);
    else
        fprintf("No solution available\n")
        return
    end
end

num_outliers_in_q = sum(z((m_normal+1):m,1));

% Calculate errors
output.totSqErrTru1 = -1;
dist = abs(X*beta_star);
if dep_var == true % get error along response direction; recall beta_1 = -1
    output.totSqErrTru1 = sum(dist(1:m_normal,1).*dist(1:m_normal,1),1);
    output.totSqErrAll1 = sum(dist(:,1).*dist(:,1),1);
    edist = dist;
    dist_sort = sort(dist(:,1).*dist(:,1));
    mio3_tsestar1 = sum(dist_sort(1:m_normal)); 
    mio3_tse1 = sum(dist_sort(1:q)); 
else % get orthogonal error
    gradLen = norm(beta_star(2:n,1)); % first coefficient is the intercept; exclude that from the gradLen calculation
    edist = abs(dist/gradLen);
    output.totSqErrTru1 = sum(edist(1:m_normal,1).*edist(1:m_normal,1));
    output.totSqErrAll1 = sum(edist(:,1).*edist(:,1));
    my_dist = edist(:,1).*edist(:,1); % squared orthogonal distances
    my_dist_sort = sort(my_dist(:,1)); % sorted squared distances
    mio3_tsestar1 = sum(my_dist_sort(1:m_normal,1));
    mio3_tse1 = sum(my_dist_sort(1:q,1));
end

if exist('errs.maxDist') > 0 
   maxDist = errs.maxDist
   output.numCloseAll1 = sum(edist <= maxDist) % count number of points within maxDist (provided on input)
end

% STEP 2: set up the PCA or OLS solution on the retained points from the MIO,
% indexed by z
%fprintf("STARTING STEP 2\n")

Xnew = zeros(m,n);
mnew = 0;
for i=1:m
   if z(i,1) % if a point is retained in the MIO solution
       mnew = mnew + 1;
       Xnew(mnew,:) = X(i,:); % add it to the new dataset
   end
end
Xnew = Xnew(1:mnew,:);


if dep_var == true % regression on the retained points
    fprintf("Using regression\n")
    beta = regress(Xnew(:,1), Xnew(:,2:n)); % regress says to keep the column of 1s to get an intercept
    beta_star = [-1.0; beta];
else % no dependent variable, so center and do PCA
    [output.w,output.w0] = getPCAHP(Xnew(:,2:n)); % omit column of ones 
    beta_star = [-output.w0; output.w(1:(n-1),1)]; % John's convention is w*x-w0=0, but Paul's is w*x+w0=0
end
beta_star = beta_star(:,1);

% Calculate errors
output.totSqErrTru2 = -1;
dist = abs(X*beta_star);
sortedabsdist = sort(dist); % for general, distance is unnormalized
gamma2 = sortedabsdist(q,1)
if dep_var == true % get error along response direction; beta_1 = -1
    output.totSqErrTru2 = sum(dist(1:m_normal,1).*dist(1:m_normal,1));
    output.totSqErrAll2 = sum(dist(:,1).*dist(:,1));
    edist = dist;
    dist_sort = sort(dist(:,1).*dist(:,1));
    mio3_tsestar2 = sum(dist_sort(1:m_normal,1)); 
    mio3_tse2 = sum(dist_sort(1:q,1)); 
else % get orthogonal error
    gradLen = norm(beta_star(2:n)); % first coefficient is the intercept; exclude that from the gradLen calculation
    edist = abs(dist/gradLen);
    output.totSqErrTru2 = sum(edist(1:m_normal,1).*edist(1:m_normal,1));
    output.totSqErrAll2 = sum(edist(:,1).*edist(:,1));
    my_dist = edist(:,1).*edist(:,1); % squared orthogonal distances
    my_dist_sort = sort(my_dist); % sorted squared distances
    mio3_tsestar2 = sum(my_dist_sort(1:m_normal,1));
    mio3_tse2 = sum(my_dist_sort(1:q,1));
end

if exist('errs.maxDist') > 0 
    output.numCloseAll2 = sum(edist <= maxDist) % count number of points within maxDist
end

%%%%%%%
% STEP 3: check which points are outliers relative to the HP found in step
% 2 and ignore those, but reinstate points that are not outliers relative
% to that HP. Find a final HP via PCA or OLS.
%fprintf("STARTING STEP 3\n")

TF = isoutlier(edist); % use distance to check which points are outliers
Xnew = zeros(m,n);
mnew = 0;
for i=1:m
   if ~TF(i,1) % if not an outlier, re-instate
       mnew = mnew + 1;
       Xnew(mnew,:) = X(i,:);
   end
end
Xnew = Xnew(1:mnew,:); % new dataset includes re-instated points

if dep_var == true
    beta = regress(Xnew(:,1), Xnew(:,2:n)); % regress says to keep the column of 1s to get an intercept
    beta_star = [-1.0; beta];
else % no dependent variable, so center and do PCA
    [output.w,output.w0] = getPCAHP(Xnew(:,2:n)); % omit column of ones 
    beta_star = [-output.w0; output.w(1:(n-1),1)]; % John's convention is w*x-w0=0, but Paul's is w*x+w0=0
end
beta_star = beta_star(:,1);

% Calculate errors
output.totSqErrTru3 = -1;
dist = abs(X*beta_star);
sortedabsdist = sort(dist); % for general, distance is unnormalized
gamma3 = sortedabsdist(q,1)
if dep_var == true % get error along response direction
    output.totSqErrTru3 = sum(dist(1:m_normal,1).*dist(1:m_normal,1));
    output.totSqErrAll3 = sum(dist(:,1).*dist(:,1));
    edist = dist;
    dist_sort = sort(dist(:,1).*dist(:,1));
    mio3_tsestar3 = sum(dist_sort(1:m_normal,1)); 
    mio3_tse3 = sum(dist_sort(1:q,1)); 
else % get orthogonal error
    gradLen = norm(beta_star(2:n)); % first coefficient is the intercept; exclude that from the gradLen calculation
    edist = abs(dist/gradLen);
    output.totSqErrTru3 = sum(edist(1:m_normal,1).*edist(1:m_normal,1));
    output.totSqErrAll3 = sum(edist(:,1).*edist(:,1));
    my_dist = edist(:,1).*edist(:,1); % squared orthogonal distances
    my_dist_sort = sort(my_dist); % sorted squared distances
    mio3_tsestar3 = sum(my_dist_sort(1:m_normal,1));
    mio3_tse3 = sum(my_dist_sort(1:q,1));
end


if exist('errs.maxDist') > 0 
    output.numCloseAll3 = sum(edist <= maxDist) % count points within maxDist
end

[out, idx] = sort(dist); % dist or edist doesn't matter because one is a constant multiple of the other
rq = edist(idx == q); % get qth distance

out_fname = strcat(resloc, "/", formulation, "i", int2str(iteration), ".csv")
disp(out_fname)
out_file = fopen(out_fname, 'w');
output.totTime  = toc;


% filename, iteration, total number of points, number of variables, number of non-outliers, percentile for LQS, formulation- mio3, 
% total squared error to hyperplane (along response or orthogonal),  gurobi runtime, gurobi status, gamma, gurobi bound,
% error for non-outliers after MIO,  error for all points after MIO, 
% error for non-outliers after PCA/OLS, error for all points after PCA/OLS, 
% error for non-outliers after step 3, error for all points after step 3, 
% num outliers in q, algorithm 3 time
% MIO1 TSE, MIO2 TSE, MIO3, TSE
% MIO1 TSEstar, MIO2 TSEstar, MIO3, TSEstar


fprintf(out_file,"%s,%d,%d,%d,%d,%d,%s,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", datafname, iteration, m, n-1, m_normal, q, formulation, output.totSqErrTru3, result.runtime, result.status, f_beta_star, result.objbound, output.totSqErrTru1, output.totSqErrAll1, output.totSqErrTru2, output.totSqErrAll2, output.totSqErrTru3, output.totSqErrAll3, num_outliers_in_q, alg3_time,mio3_tse1,mio3_tse2,mio3_tse3,mio3_tsestar1,mio3_tsestar2,mio3_tsestar3,gamma2,gamma3,output.totTime);

fclose(out_file);

return
end
