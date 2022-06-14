% May 26, 2022
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% Finds a best fitting regressionhyperplane using mixed-integer optimization 
% and an advanced start provided by CBreg. The steps are:
% 1. Run CBreg to obtain heuristic solution
% 2. MIO1: run the improved MIO minimization of the qth percentile error, 
%    in the sorted set of errors, subject to the time limit.
% 3. MIO2: place a PCA hyperplane on the q points kept by the MIO solution.
% 4. MIO3: check which points are outliers relative to the step 2 HP. Outliers
%    removed, nonoutliers kept, and a final HP is placed. Generally there
%    is likely to be only reintroduction of points that are not outliers.
%
%NOTE: The regression equation is w0 + w1x1 + w2x2 + ... + wnxn = y,
% so a column of ones is added at the beginning of the input Ain matrix
%
%INPUTS:
%  y: the output variable values
%  Ain: the data matrix (points x predictor variables)
%  mioparams: struct containing the CB and MIO control parameters:
%    .mtru: the number of inliers, if known. Useful in generating
%      statistics for testing purposes. Not used if < 1. If used, then the
%      first mtru points in the data set must be the inliers.
%    .q: MIO minimizes the qth deviation from the hyperplane. 
%      - if q > 0 it is the percentile of the m points in the data set
%      - if q < 0 then -q is the actual ordinal number (not percentile)
%      - if q = 0 then take the q from CBgen
%    .maxResid: points with residual error less than this are "close". 
%      There are 3 cases:
%         < 0: -maxResid is the percentile of residuals for the first
%               hyperplane. Note it is in %. 
%               NOTE: maxResid = -16 IS HIGHLY RECOMMENDED.
%               If maxResid is not specified, then -16 is used.
%         = 0: means that closeAll is not used to help identify the
%              best hyperplane. Result is just the final hyperplane.
%         > 0: an actual residual to define maxResid. 
%  gbparams: struct containing all of the Gurobi control parameters
%    Sample values:
%      .TimeLimit = 3600;
%      .OutputFlag = 1;
% OUTPUTS:
%  result: the Gurobi result struct
%  output: the fitted hyperplane output struct, with these fields:
%    .w: the weights for the fitted hyperplane
%    .gamma(k,1): the minimum value of the residual at the qth percentile
%    .TSEstar(k,1): if mtru > 0, the sum of the squared residuals
%       for the mtru smallest residuals, for each of the steps, plus output:
%         k=1: CBgen result, if used
%         k=2: MIO result
%         k=3: PCA on MIO output
%         k=4: reinstatement, then final PCA
%    .TSE(k,1): the sum of the squared residuals for the q points 
%       having the smallest residuals, for each of the steps, plus output
%    IF maxResid ~= 0, then these values are also calculated:
%      .closeAll(k,1): the number of points having a residual of 
%        less than maxResid
%      .bnd2: sum of squared residuals to a regression hyperplane fit to the
%        mtru input points
%
% DEPENDENCIES: this routine calls CBreg, and the external solver Gurobi.
%  Gurobi has a free academic license.

function [result,output] = CBMIOreg(y,Ain,mioparams,gbparams)

mtru = mioparams.mtru;

% if unused, values are -1.
output.TSEstar = zeros(4,1) - 1;
output.TSE = zeros(4,1) - 1;
output.closeAll = zeros(4,1) - 1;
output.gamma = zeros(4,1) - 1;
output.maxResid = -1;

tStart = tic;

% The regression equation is w0 + w1x1 + w2x2 + ... + wnxn = y,
% so a column of ones is added at the beginning of Ain for w0
A = [ones(size(Ain,1),1),Ain];
m = size(A,1);
n = size(A,2);

% set an integer value of q
q = mioparams.q;
if q > 0
    % input q is a percentile
    q = floor(q/100*m);
else
    if q < 0
        % input q is a specific integer value
        q = -q;
    end
end

%--------------------------------------------------------------------------
% STEP 1: run CBreg

fprintf("-----CBreg starts-----\n")
CBparam.maxResid = mioparams.maxResid;
CBparam.mgood = mioparams.mtru;
[inc] = CBreg(y,Ain,CBparam);
fprintf("-----CBreg ends-----\n")
maxResid = inc.maxResid;
output.maxResid = maxResid;
output.bnd2 = inc.bnd2;
fprintf("  CBreg sets maxResid at %f\n",maxResid)

if q == 0
    % Set q based on outFinder results, found in CBreg ouput
    q = inc.q;
end

output.q = q;
fprintf("  q set at %d\n",q)
[output] = update(1,y,Ain,inc.w0Out,inc.weightsOut,q,maxResid,mtru,output);
output.CBregTime = toc(tStart);

if inc.status == -1
    fprintf("  CBreg failure: aborting MIO solution.\n")
    result.status = "ABORTED";
    result.mipgap = Inf;
    return
end
if inc.status == 1
    fprintf("  CBreg exact solution: aborting MIO solution.\n")
    result.status = "ABORTED";
    result.mipgap = Inf;
    return
end

%--------------------------------------------------------------------------
% STEP 2: set up and solve the MIO 
fprintf("Solving the MIO to minimize gamma over all %d points.\n",m)
tMIOstart = tic;

% Initialize the constraint matrix
% row order: m elastic cons, m error cons, 1 summary q constraint
% col order: n w, m eplus, m eminus, m rel, m z, 1 gamma
model.A = sparse(2*m+1,n+4*m+1);
model.rhs = zeros(2*m+1,1);
% m Elastic cons, of form w0 + w1x1 + w2x2 + ... +wnxn + eplus - eminus = y. 
output.w = zeros(n,1);
model.A(1:m,1:n+2*m) = [sparse(A),speye(m),-speye(m)];
model.rhs(1:m,1) = y;
model.sense(1:m,1) = repmat('=',m,1);
% m error constraints, of form eplus_i + eminus_i - rel_i - gamma < 0
model.A(m+1:2*m,n+1:n+3*m) = [speye(m),speye(m),-speye(m)];
model.A(m+1:2*m,n+4*m+1) = sparse(-ones(m,1));
model.sense(m+1:2*m,1) = repmat('<',m,1);
% the q equation, of form sum(Z) = q
model.A(2*m+1,n+3*m+1:n+4*m) = sparse(ones(1,m));
model.rhs(2*m+1,1) = q;
model.sense(2*m+1,1) = '=';
% variable lower bounds
model.lb = [zeros(n,1)-Inf;zeros(m*4+1,1)];
% variable upper bounds
model.ub = zeros(n+4*m+1,1)+Inf;
model.ub(n+3*m+1:n+4*m,1) = ones(m,1);
% set up the SOS constraints
for i=1:m
    % rel,z
    model.sos(i).type = 1;
    model.sos(i).index = [n+2*m+i,n+3*m+i];
end
% specify variable types
model.vtype = [repmat('C',n+3*m,1);repmat('B',m,1);'C'];
% set up objective function to minimize gamma
model.obj = [zeros(n+4*m,1);1.0];
model.modelsense = 'min';

% Use the CBreg advanced start
% Initialize some variables
eplus = zeros(m,1);
eminus = zeros(m,1);
rel = zeros(m,1);
z = zeros(m,1);
% Calculate gamma for CBreg eqn
sortedabsresid = sort(abs(inc.residOut));
gamma = sortedabsresid(q,1);

for i=1:m
   if inc.residOut(i,1) > 0
       eminus(i,1) = inc.residOut(i,1);
   else
       eplus(i,1) = -inc.residOut(i,1);
   end
   if abs(inc.residOut(i,1)) <= gamma
       z(i,1) = 1;
       rel(i,1) = 0;
   else
       rel(i,1) = abs(inc.residOut(i,1));
   end
end

model.StartNumber = 0;
model.start = [inc.w0Out;inc.weightsOut; eplus; eminus; rel; z; gamma];

% solve the model
result = gurobi(model, gbparams);

% About solutions:
% Check results.status:
%   If 'OPTIMAL' then take result.x
%   If 'TIME_LIMIT' then check whether result.mipgap is Inf
%     If yes, exit with no solution.
%     Else take result.pool.objval(1,1) as incumbent value of objective and
%       take result.pool(1).xn as the associated solution.

% get results
fprintf("  MIO solution status: %s\n",result.status)
if strcmp(result.status, 'OPTIMAL')
    mioOut.w0 = result.x(1,1);
    mioOut.w = result.x(2:n,1);
%     eplus = result.x(n+1:n+m,1);
%     eminus = result.x(n+m+1:n+2*m,1);
%     rel = result.x(n+2*m+1:n+3*m,1);
    z = result.x(n+3*m+1:n+4*m,1);
    mioOut.gamma = result.x(n+4*m+1,1);
else
    if result.mipgap ~= Inf
        % Gurobi stopped for some reason, maybe time limit, but it has an
        % incumbent solution, so use that
        fprintf("  Using incumbent solution\n")
        mioOut.w0 = result.pool(1).xn(1,1);
        mioOut.w = result.pool(1).xn(2:n,1);
%         eplus = result.pool(1).xn(n+1:n+m,1);
%         eminus = result.pool(1).xn(n+m+1:n+2*m,1);
%         rel = result.pool(1).xn(n+2*m+1:n+3*m,1);
        z = result.pool(1).xn(n+3*m+1:n+4*m,1);
        mioOut.gamma = result.pool(1).xn(n+4*m+1,1);
    else
        fprintf("  No MIO solution available: aborting.\n")
        result.status = strcat(result.status,'_and_MIO_failure');
        return
    end
end
fprintf("  MIO gamma %f at q = %d\n",mioOut.gamma,q)

% Update solution
[output] = update(2,y,Ain,mioOut.w0,mioOut.w,q,maxResid,mtru,output);
output.MIOTime = toc(tMIOstart);
fprintf("  MIO solution time %f\n",output.MIOTime)

%-------------------------------------------------------------------------
% STEP 3: multiple regression solution on the retained points from the MIO,
% indexed by z

Anew = zeros(m,n-1);
ynew = zeros(m,1);
mnew = 0;
for i=1:m
   if z(i,1)
       mnew = mnew + 1;
       Anew(mnew,:) = Ain(i,:);
       ynew(mnew,1) = y(i,1);
   end
end
Anew = Anew(1:mnew,:);
ynew = ynew(1:mnew,1);
Anew = [ones(mnew,1),Anew];
fprintf("Running multiple regression on %d points retained after MIO.\n",mnew)
beta = regress(ynew,Anew);
w0 = beta(1,1);
w = beta(2:n,1);

% Update solution
[output] = update(3,y,Ain,w0,w,q,maxResid,mtru,output);

%--------------------------------------------------------------------------
% STEP 4: check which points are outliers relative to the HP found in step
% 3 and ignore those, but reinstate points that are not outliers relative
% to that HP. Find a final HP via multiple regression.

% Calculate residuals
resid = w0 + Ain*w - y;
absresid = abs(resid);
TF = isoutlier(absresid);
Anew = zeros(m,n-1);
ynew = zeros(m,1);
mnew = 0;
for i=1:m
   if ~TF(i,1)
       mnew = mnew + 1;
       Anew(mnew,:) = Ain(i,:);
       ynew(mnew,1) = y(i,1);
   end
end
Anew = Anew(1:mnew,:);
Anew =[ones(mnew,1),Anew];
ynew = ynew(1:mnew,1);
fprintf("Reinstating/removing points relative to last hyperplane, and rerunning multiple regression on %d pts\n",mnew)
beta = regress(ynew,Anew);
w0 = beta(1,1);
w = beta(2:n,1);

% Update solution
[output] = update(4,y,Ain,w0,w,q,maxResid,mtru,output);

output.w0 = w0;
output.w = w;
output.solTime  = toc(tStart);
fprintf("Solution time: %f\n",output.solTime);

return
end

%--------------------------------------------------------------------------
% Make calculations and updates at step k
function [output] = update(k,y,A,w0,w,q,maxResid,mtru,output)
% Calculate the residuals
resid = w0 + A*w - y;
absresid = abs(resid);
sortedresid = sort(absresid);
output.gamma(k,1) = sortedresid(q,1);
fprintf("  gamma %f at q = %d\n",output.gamma(k,1),q)
% Calculate TSE
output.TSE(k,1) = norm(sortedresid(1:q,1).*sortedresid(1:q,1),1);
fprintf("  TSE %f at q = %d\n",output.TSE(k,1),q);
if mtru > 0
    % Calculate TSEstar
    output.TSEstar(k,1) = norm(sortedresid(1:mtru,1).*sortedresid(1:mtru,1),1);
    fprintf("  TSEstar %f at mtru = %d\n",output.TSEstar(k,1),mtru);
end
if maxResid > 0
    output.closeAll(k,1) = sum(absresid <= maxResid);
    fprintf("  %d close points\n",output.closeAll(k,1));
    output.TSEstar(k,1) = norm(sortedresid(1:mtru,1).*sortedresid(1:mtru,1),1);
end

return
end
