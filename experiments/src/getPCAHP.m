%==========================================================================
% February 27, 2024
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA
%
% Get the best-fitting hyperplane from an input dataset, using PCA.
% If there are insufficient points, default to an elastic linear
% programming solution, which uses the Gurobi LP solver.
%
% INPUTS: dataSet is the data matrix
% OUTPUTS: Hyperplane equation is: w*x = RHS.
%   w: weights in the hyperplane equation
%   RHS: right hand side constant in the hyperplane equation
function [w,RHS] = getPCAHP(dataSet)

n = size(dataSet,2);

% pca doesn't work if the number of rows in the dataset is less than the
% number of columns, so run an elastic LP solution instead
if size(dataSet,1) < n
    fprintf("  Too few observations vs. features for PCA: running elastic LP.\n")
    m = size(dataSet,1);
    % Set Gurobi parameters
    gbparams = struct();
    gbparams.OutputFlag = 0;
    gbparams.Threads=14;
    % Build the complete model constraint matrix and solve.
    model.A = sparse([dataSet,speye(m,m),-speye(m,m)]);
    model.rhs = zeros(m,1) + n;
    model.sense(1:m,1) = repmat('=',m,1);
    nfinal = n + 2*m;
    % variable lower bounds
    model.lb = zeros(nfinal,1);
    model.lb(1:n,1) = -Inf;
    % variable upper bounds
    model.ub = zeros(nfinal,1) + Inf;
    % specify variable types
    model.vtype = repmat('C',nfinal,1);
    % set up objective function
    model.obj = ones(nfinal,1);
    model.obj(1:n,1) = 0.0;
    model.modelsense = 'min';
    % Solve LP to find hyperplane
    result = gurobi(model, gbparams);
    w = result.x(1:n,1);
    RHS = n;
    return
end

% column means
my_means = mean(dataSet);
% pca on centered data
my_loadings = pca(dataSet,'Economy',false);
% intercept of hyperplane
RHS = sum(my_means*my_loadings(:,n));
% norm vector of hyperplane
w = my_loadings(:,n);
% Hyperplane equation is: w*x = RHS.
return

end
