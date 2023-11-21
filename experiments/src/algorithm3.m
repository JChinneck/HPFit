% Algorithm 3 is proposed in Bertsimas and Mazumder (2014).  The 
% result is used as a warm start for a mixed-integer optimization
% approach.  
% algorithm3() is a composite heuristic for fitting a hyperplane to 
% minimize the q^th residual.  This version is called by mio.m. 
% To run algorithm 3 by itself, use run_alg3.m
% Purpose: to find the best fitting hyperplance to a set of input 
% data, where best fitting means that the absolute value of the 
% q^th residual is minimized.  Bertsimas and Mazumder (2014) call
% this Least Quantile of Squres (LQS) even though there is no square.
% which can be affected by outliers.
% Algorithm 3 begins by running L1 regression (LAD).  If no response
% variable is specified, we use PCA. 
% The result is
% perturbed 100 times and sent to Algorithm 2 for 500 iterations. 
% Algorithm 2 is a subdifferential-based algorithm for LQS.  The best
% solution is passed to Algorithm 1.  Algorithm 1 is a sequential 
% linear optimization algorithm for LQS.
% Algorithms 1 and 2 have been adapted to the case when no response
% variable is specified. 

% Main options: 
% q: 
%   - the order statistic used to evaluate the fit of a hyperplane. 
%     The fit is determined by the q^th largest residual.
% dep_var:
%   - true: a dependent variable is specified, as in ordinary
%           regression.  A residual is measured as the absolute
%           difference between the response value and the 
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent variable is specified. The intercept term 
%            is fixed to n, the original number of variables in the
%            dataset.  A residual is measured using the elastic LP
%            measure which is the absolute difference between 
%            beta^T x (including beta_0=n) and zero.
% 
% NOTES:
% - when a dependent variable is specified, the corresponding 
%   coefficient is set to -1.  When one is not specified, the 
%   intercept is set to n.
%
% INPUTS: 
% - options mentioned above
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

function [beta_star, f_beta_star] = algorithm3(X,q,dep_var,init_method) % starts near the bottom; subfunctions defined at the top
rng('default')
%Algorithm 3 starts here -------------------------------------------------
%X  = X{:,:};
[m,n] = size(X); % get the size of dataset; for datasets with no response, a column of 1s is included
beta = lad(X, dep_var, init_method); % get LAD solution
f_beta_star = inf;
nu = 2;

for k=1:100 % perturb LAD solution 100 times
    if dep_var == true
        perturb = -nu + 2*nu*rand(n-1,1); % random number in [-nu, nu]
    else
        perturb = -nu + 2*nu*rand(n,1); % random number in [-nu, nu]
    end
    beta = beta + times(perturb, abs(beta));
    [f_beta, beta] = algorithm2(X, beta, 500, q, dep_var); % use perturbed LAD as starting point for Algorithm 2
    if f_beta < f_beta_star % if better, update
        beta_star = beta;
        f_beta_star = f_beta;
    end
end
[beta_star, f_beta_star] = algorithm1(X, beta_star, q, dep_var); % use best Algorithm 2 solution as input to Algorithm 1

if dep_var == true
    beta_star = [-1;beta_star]; % with a dependent variable, add the coefficient back on y
else
    beta_star
    beta_star = beta_star*((n-1)/beta_star(1,1))
end


end

function [beta_star] = lad(X, dep_var, init_method) % LAD = least aboslute deviations = L1 regression
    [m,n] = size(X); 
    f_beta_star = inf;

    if dep_var == true % dependent variable - L1 regression
       model.obj = [zeros(n-1,1) ; ones(2*m,1) ];  % beta; eplus; eminus
       model.lb  = [-inf(n-1,1) ; zeros(2*m,1) ];
       model.A   = [ sparse(X(:,2:n)) speye(m) -speye(m)]; % beta^Tx + eplus - eminus = 0 for each point i
       model.rhs = [ X(:,1) ]; % RHS is response variable
       model.sense = repmat('=', m, 1);
       model.modelsense = 'min';
       params = struct();
       %params.OutputFlag = 0;
       result = gurobi(model, params);
       beta_star = result.x(1:(n-1)); % coefficients for variables except for response, which is -1 (we don't store the response coefficient)
       f_beta_star = result.objval;
    else % no dependent variable
        if init_method == "LP"  % elastic LP
            model.obj = [zeros(n-1,1) ; ones(2*m,1) ]; % beta (intercept fixed at n), eplus, eminus
            model.lb  = [-inf(n-1,1) ; zeros(2*m,1) ];
            model.A   = [ sparse(X(:,2:n)) speye(m) -speye(m)]; % beta^T X + eplus - eminus = n, for each point i
            model.rhs = [ zeros(m,1) - n ];
            model.sense = repmat('=', m, 1);
            model.modelsense = 'min';
            params = struct();
            params.OutputFlag = 0;
            result = gurobi(model, params);
            beta_star = [n;result.x(1:n-1)];
        else % PCA
            [weights, intercept] = getPCAHP(X(:,2:n)); % use PCA to get weights and intercept
            beta_star = [-intercept; weights]; % put the intercept first to match up with the column of 1s
        end
    end
end

function [beta_star] =solve_lp(X, beta, q, dep_var) % called by Algorithm 1
    [m,n] = size(X); 
    if dep_var == true % first variable is response
        residuals = X(:,1) - X(:,2:n)*beta; % residuals are responses minus projections on hyperplane
        abs_residuals = abs(residuals);
        [out,idx] = sort(abs_residuals); % sort absolute residuals and get indexes
        w_star = zeros(m,1);
        w_star(idx(m-q+1:m),1) = 1; % get q largest residuals
        subg = transpose(X(:,2:n))*(times(-w_star,sign(residuals))); % calculate subgradient
        
        model.obj = [ (m-q+1); ones(m, 1); -subg ]; % theta; nu; beta
        model.lb  = [ -inf ; zeros(m,1); -inf(n-1,1) ];
        model.A   = [ sparse(ones(m, 1)) speye(m) sparse(X(:,2:n));  % theta + nu_i >= y_i - x_i^T beta, for each point i
            sparse(ones(m,1)) speye(m) sparse(-X(:,2:n)); % theta + nu_i >= -y_i+x_i^T beta for each point i
            (m-q+1) sparse(ones(1,m)) -subg' ]; % objective is at least 0
        model.rhs = [X(:,1); -X(:,1); 0.0];  % +/- y_i
        model.sense = [repmat('>', 2*m, 1); '>'];
        model.modelsense = 'min';
        params = struct();
        params.OutputFlag = 0;
        result = gurobi(model, params);
        beta_star = result.x((m+2):(m+1+n-1)); % beta
    else % no dependent variable
        residuals = X*beta;
        abs_residuals = abs(residuals);
        [out,idx] = sort(abs_residuals);
        w_star = zeros(m,1);
        w_star(idx(m-q+1:m),1) = 1; % get q largest residuals
        subg = transpose(X)*(times(-w_star,sign(residuals)));
        
        model.obj = [ (m-q+1); ones(m, 1); -subg ]; % theta; nu_i; beta
        model.lb  = [ -inf ; zeros(m,1); -inf(n,1) ];
        model.A   = [ sparse(ones(m, 1)) speye(m) sparse(X); % theta + nu_i >= -x^T beta 
        sparse(ones(m,1)) speye(m) sparse(-X); % theta + nu_i >= x^T beta
        sparse(zeros(1,1)) sparse(zeros(1,m)) 1.0 sparse(zeros(1,n-1)); % beta_0 = n (n-1 because of column of ones for constant)
        (m-q+1) sparse(ones(1,m)) -subg' ]; % objective is at least 0
        model.rhs = [zeros(m,1); zeros(m,1) ; n-1; 0.0]; % 0; 0; n-1; 0
        model.sense = [repmat('>', 2*m, 1) ; '='; '>'];
        model.modelsense = 'min';
        params = struct();
        params.OutputFlag = 0;
        result = gurobi(model, params);
        beta_star = result.x((m+2):(m+1+n));
    end
end

function [beta_kplus1, q_residual_k] = algorithm1(X, beta, q, dep_var) % sequential linear optimization algorithm.  Last part of Algorithm 3.
    [m,n] = size(X); 
    Tol = 0.0001;
    k=0;
    if dep_var == true
        residuals_k = abs(X(:,1) - X(:,2:n)*beta); % absolute value of response minus projection on hyperplane
        [out_k, idx_k] = sort(residuals_k);
        q_residual_k = residuals_k(idx_k(q)); % q^th largest absolute residual
        while 0 < 1
            beta_kplus1 = solve_lp(X, beta, q, dep_var); % get new beta
            residuals_kplus1 = abs(X(:,1) - X(:,2:n)*beta_kplus1); % calculate new residuals
            [out_kplus1, idx_kplus1] = sort(residuals_kplus1); % sort residuals and get indexes
            q_residual_kplus1 = residuals_kplus1(idx_kplus1(q)); % q^th largest residual
            if abs(q_residual_k - q_residual_kplus1) <= Tol*q_residual_k % if converged, stop
                return
            end
            q_residual_k = q_residual_kplus1;
            k = k+1;
        end
    else % no dependent variable
        residuals_k = abs(X*beta); % absolute value of distance to hyperplane using elastic LP criterion
        [out_k, idx_k] = sort(residuals_k); % sort residuals and get indexes
        q_residual_k = residuals_k(idx_k(q)); % q^th largest residual
        while 0 < 1
            beta_kplus1 = solve_lp(X, beta, q, dep_var); % get new beta
            residuals_kplus1 = abs(X*beta_kplus1); % calculate new residuals
            [out_kplus1, idx_kplus1] = sort(residuals_kplus1); % sort new residuals
            q_residual_kplus1 = residuals_kplus1(idx_kplus1(q)); % q^th largest residual
            if abs(q_residual_k - q_residual_kplus1) <= Tol*q_residual_k % if converged, stop
                return
            end
            q_residual_k = q_residual_kplus1;
            k = k+1;
        end
    end
end 

function [f_beta_star, beta_star] = algorithm2(X, beta, MaxIter, q, dep_var) % subgradient optimization algorithm
    beta_star = beta;
    [m,n] = size(X);
    if dep_var == true % if a response variable is specified as in regression
        l2_norms = vecnorm(X(:,2:n),2,2); % get 2-norm of each row
        alpha_k = 1/max(l2_norms); % stepsize is constant as in paper
        residuals = X(:,1) - X(:,2:n)*beta_star; % residuals are response minus projections on hyperplane
        abs_residuals = abs(residuals);
        [out, idx] = sort(abs_residuals); % sort absolute residuals and get indexes
        f_beta_star = abs_residuals(idx(q)); % get q^th largest residual
        for k = 1:MaxIter
            beta_kplus1 = beta - alpha_k * (-sign(residuals(idx(q))))*transpose(X(idx(q),2:n)); % update beta
            residuals = X(:,1) - X(:,2:n)*beta_kplus1; % calculate new residuals
            abs_residuals = abs(residuals);
            [out, idx] = sort(abs_residuals); % sort residuals and get indexes
            f_beta = abs_residuals(idx(q)); % get q^th largest residual
            sortedabsresid = sort(abs_residuals);
            if f_beta < f_beta_star % if the q^th largest is better, update the solution
                f_beta_star = f_beta;
                beta_star = beta_kplus1;
            end
            beta = beta_kplus1;
        end
    else % no dependent variable
        l2_norms = vecnorm(X,2,2); % get 2-norm of each row
        alpha_k = 1/max(l2_norms); % stepsize is constant as in paper
        residuals = X*beta_star; % calculate residuals using elastic LP criterion
        abs_residuals = abs(residuals);
        [out, idx] = sort(abs_residuals); % sort residuals and get indexes
        f_beta_star = abs_residuals(idx(q)); % get q^th largest residual
        for k = 1:MaxIter
            beta_kplus1 = beta - alpha_k * (-sign(residuals(idx(q))))*transpose(X(idx(q),:)); % update beta
            residuals = X*beta_kplus1; % calculate new residuals
            abs_residuals = abs(residuals);
            [out, idx] = sort(abs_residuals); % sort new residuals and get indexes
            f_beta = abs_residuals(idx(q)); % get q^th largest residual
            if f_beta < f_beta_star % if q^th largest is better, update
                f_beta_star = f_beta;
                beta_star = beta_kplus1;
            end
            beta = beta_kplus1;
        end
    end
    return
end
