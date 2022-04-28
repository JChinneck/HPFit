% Get the best-fitting hyperplane from an input dataset, using PCA
function [w,RHS] = getPCAHP(dataSet)

% Code from Paul Brooks
n = size(dataSet,2);

% pca doesn't work if the number of rows in the dataset is less than the
% number of columns, so run an elastic LP solution instead
if size(dataSet,1) < n
    fprintf("  too few observations vs. features for pca to run: running elastic LP\n")
    m = size(dataSet,1);
    % Set MOSEK parameters
    param.MSK_IPAR_LOG = 0;
    param.MSK_IPAR_OPF_WRITE_HEADER = 'MSK_OFF';
    param.MSK_IPAR_WRITE_SOL_HEAD = 'MSK_OFF';
    % Build the complete model constraint matrix and solve.
    Alp = sparse([dataSet,speye(m,m),-speye(m,m)]);
    blc = zeros(m,1) + n;
    buc = zeros(m,1) + n;
    nfinal = size(Alp,2);
    blx = zeros(nfinal,1);
    bux = zeros(nfinal,1) + Inf;
    blx(1:n,1) = -Inf;
    % set up original objective function
    c = ones(nfinal,1);
    c(1:n,1) = 0.0;
    c(n+1:n+m,1) = ones(m,1);
    c(n+m+1:n+2*m,1) = ones(m,1);
    % Solve LP to find hyperplane
    [res] = msklpopt(c,Alp,blc,buc,blx,bux,param,'minimize echo(0)');
    w = res.sol.bas.xx(1:n,1);
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
% eqn is: pca_norm*x - pca_intercept = 0.
return

end
