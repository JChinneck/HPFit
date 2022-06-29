% June 28, 2022
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% Finds a best fitting hyperplane using mixed-integer optimization and an 
% advanced start provided by CBgen. The steps are:
% 1. Run CBgen to obtain heuristic solution
% 2. Run the improved MIO minimization of the qth percentile error, in the
%    sorted set of errors, subject to the time limit.
% 3. Places a PCA hyperplane on the q points kept by the MIO solution.
% 4. Checks which points are outliers relative to the step 2 HP. Outliers
%    removed, nonoutliers kept, and a final HP is placed. Generally there
%    is likely to be only reintroduction of points that are not outliers.
%
%NOTE: hyperplane equations have the form: w_1*x_1 + w_2*x_2 + .... = RHS
%
%INPUTS:
%  Ain: the data matrix
%  mioparams: struct containing the CB and MIO control parameters:
%    .mtru: the number of inliers, if known. Useful in generating
%      statistics for testing purposes. Not used if < 1. If used, then the
%      first mtru points in the data set must be the inliers.
%    .q: MIO minimizes the qth deviation from the hyperplane. 
%      - if q > 0 it is the percentile of the m points in the data set
%      - if q < 0 then -q is the actual ordinal number (not percentile)
%      - if q = 0 then take the q from CBgen
%    .maxDist: points closer than this distance to a hyperplane are
%      "close". Distance is Euclidean.
%      There are 3 cases:
%         < 0: -maxDist is the percentile of distances from the first
%               hyperplane. Note it is in %. 
%               NOTE: maxDist = -16 IS HIGHLY RECOMMENDED.
%               If maxDist is not specified, then -16 is used.
%         = 0: means that closeAll is not used to help identify the
%              best hyperplane. Result is just the final hyperplane.
%         > 0: an actual Euclidean distance to define maxDist. 
%  gbparams: struct containing all of the Gurobi control parameters
%    Sample values:
%      .TimeLimit = 3600;
%      .OutputFlag = 1;
% OUTPUTS:
%  result: the Gurobi result struct
%  output: the fitted hyperplane output struct, with these fields:
%    .w: the weights for the fitted hyperplane
%    .RHS: value of constant in the hyperplane expression
%    .gamma(k,1): the minimum value of the Euclidean error at the qth percentile
%    .TSEstar(k,1): if mtru > 0, the sum of the squared Euclidean distances
%       for the mtru points closest to the fitted hyperplane, for each of 
%       the steps, plus output:
%         k=1: CBgen result, if used
%         k=2: MIO result
%         k=3: PCA on MIO output
%         k=4: reinstatement, then final PCA
%    .TSE(k,1): the sum of the squared Euclidean distances for the q points 
%       closest to the fitted hyperplane, for each of the steps, plus output
%    IF maxDist ~= 0, then these values are also calculated:
%      .closeAll(k,1): the number of points within a Euclidean distance of 
%        less than maxDist from the hyperplane
%      .bnd2: sum of squared Euclidean distances to a hyperplane fit to the
%        mtru input points
%
% DEPENDENCIES: this routine calls CBgen, and the external solver Gurobi.
%    Gurobi has a free academic license.

function [result,output] = CBMIOgen(Ain,mioparams,gbparams)

mtru = mioparams.mtru;

% if unused, values are -1.
output.TSEstar = zeros(4,1) - 1;
output.TSE = zeros(4,1) - 1;
output.closeAll = zeros(4,1) - 1;
output.gamma = zeros(4,1) - 1;
output.maxDist = -1;

tStart = tic;

% get m, n, and specific q
m = size(Ain,1);
n = size(Ain,2);

% set an integer value of q
q = mioparams.q;
if isa(q, 'double')
    if q > 0
        % input q is a percentile
        q = floor(q/100*m);
    else
        if q < 0
            % input q is a specific integer value
            q = -q;
        end
    end
end

%--------------------------------------------------------------------------
% STEP 1: run CBgen

fprintf("-----CBgen starts-----\n")
CBparam.maxDist = mioparams.maxDist;
CBparam.mgood = mioparams.mtru;
[inc] = CBgen(Ain,CBparam);
fprintf("-----CBgen ends-----\n")
maxDist = inc.maxDist;
output.maxDist = maxDist;
output.bnd2 = inc.bnd2;
fprintf("  CBgen sets maxDist at %f\n",maxDist)

if isa(q, 'double')
    if q == 0
        % Set q based on outFinder results, found in CBgen ouput
        q = inc.q;
    end
end
if isa(q, 'char')
    if q == "qout"
        q = inc.qout;
    end
end
output.q = q;
fprintf("  q set at %d\n",q)

output.cbq = inc.q;
output.outfinderq = inc.qout;

[output] = update(1,Ain,inc.weightsOut,inc.RHSOut,q,maxDist,mtru,output);
output.CBgenTime = toc(tStart);

if inc.status == -1
    fprintf("  CBgen failure: aborting MIO solution.\n")
    result.status = "ABORTED";
    result.mipgap = Inf;
    return
end
if inc.status == 1
    fprintf("  CBgen exact solution: aborting MIO solution.\n")
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
% m elastic cons, of form w1x1 + w2x2 + ... +wnxn + eplus - eminus = n. 
% Note use of the RHS output by CBgen
output.w = zeros(n,1);
output.RHS = inc.RHSOut;
model.A(1:m,1:n+2*m) = [sparse(Ain),speye(m),-speye(m)];
model.rhs(1:m,1) = zeros(m,1) + inc.RHSOut;
model.sense(1:m,1) = repmat('=',m,1);
% m error constraints of form eplus_i + eminus_i - rel_i - gamma < 0
model.A(m+1:2*m,n+1:n+3*m) = [speye(m),speye(m),-speye(m)];
model.A(m+1:2*m,n+4*m+1) = sparse(-ones(m,1));
model.sense(m+1:2*m,1) = repmat('<',m,1);
% the q equation, of form sum(z) = q
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

% Use the CBgen advanced start
% Initialize some variables
eplus = zeros(m,1);
eminus = zeros(m,1);
rel = zeros(m,1);
z = zeros(m,1);
% Calculate dists with new HP eqn
dist = inc.RHSOut - Ain*inc.weightsOut;
absdist = abs(dist);
sortedabsdist = sort(absdist);
gamma = sortedabsdist(q,1);
for i=1:m
   if dist(i,1) > 0
       eplus(i,1) = dist(i,1);
   else
       eminus(i,1) = -dist(i,1);
   end
   if absdist(i,1) <= gamma
       z(i,1) = 1;
       rel(i,1) = 0;
   else
       rel(i,1) = absdist(i,1);
   end
end
model.StartNumber = 0;
model.start = [inc.weightsOut; eplus; eminus; rel; z; gamma];

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
    mioOut.w = result.x(1:n,1);
    mioOut.RHS = inc.RHSOut;
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
        mioOut.w = result.pool(1).xn(1:n,1);
        mioOut.RHS = inc.RHSOut;
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
[output] = update(2,Ain,mioOut.w,mioOut.RHS,q,maxDist,mtru,output);
output.MIOTime = toc(tMIOstart);
fprintf("  MIO solution time %f\n",output.MIOTime)

%-------------------------------------------------------------------------
% STEP 3: PCA solution on the retained points from the MIO, indexed by z

Anew = zeros(m,n);
mnew = 0;
for i=1:m
   if z(i,1)
       mnew = mnew + 1;
       Anew(mnew,:) = Ain(i,:);
   end
end
Anew = Anew(1:mnew,:);
fprintf("Running PCA on %d points retained after MIO.\n",mnew)
[weights,RHS] = getPCAHP(Anew);

% Update solution
[output] = update(3,Ain,weights,RHS,q,maxDist,mtru,output);

%--------------------------------------------------------------------------
% STEP 4: check which points are outliers relative to the HP found in step
% 3 and ignore those, but reinstate points that are not outliers relative
% to that HP. Find a final HP via PCA.

% Calculate edist
dist = Ain*weights - RHS;
gradLen = norm(weights);
edist = abs(dist/gradLen);
TF = isoutlier(edist);
Anew = zeros(m,n);
mnew = 0;
for i=1:m
   if ~TF(i,1)
       mnew = mnew + 1;
       Anew(mnew,:) = Ain(i,:);
   end
end
Anew = Anew(1:mnew,:);
fprintf("Reinstating/removing points relative to last hyperplane, and rerunning PCA on %d pts\n",mnew)

[weights,RHS] = getPCAHP(Anew);

% Update solution
[output] = update(4,Ain,weights,RHS,q,maxDist,mtru,output);
output.w = weights;
output.RHS = RHS;
output.solTime  = toc(tStart);
fprintf("Solution time: %f\n",output.solTime);

return
end

%--------------------------------------------------------------------------
% Make calculations and updates at step k
function [output] = update(k,A,weights,RHS,q,maxDist,mtru,output)
% Calculate the point distances from the hyperplane
dist = A*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);
sortededist = sort(edist);
output.gamma(k,1) = sortededist(q,1);
fprintf("  gamma %f at q = %d\n",output.gamma(k,1),q)
% Calculate TSE
output.TSE(k,1) = norm(sortededist(1:q,1).*sortededist(1:q,1),1);
fprintf("  TSE %f at q = %d\n",output.TSE(k,1),q);
if mtru > 0
    % Calculate TSEstar
    output.TSEstar(k,1) = norm(sortededist(1:mtru,1).*sortededist(1:mtru,1),1);
    fprintf("  TSEstar %f at mtru = %d\n",output.TSEstar(k,1),mtru);
end
if maxDist > 0
    output.closeAll(k,1) = sum(edist <= maxDist);
    fprintf("  %d close points\n",output.closeAll(k,1));
    output.TSEstar(k,1) = norm(sortededist(1:mtru,1).*sortededist(1:mtru,1),1);
end

return
end
