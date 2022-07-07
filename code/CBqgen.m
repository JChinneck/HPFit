% July 6, 2022
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% Heuristic function to minimize the qth error in a general data set. Main
% steps are:
%   1. Run CBgen to obtain hyperplane fit.
%   2. Select the q points closest to the the CBgen HP
%   3. Solve an LP on the q retained points to minimize the maximum error
%   4. Return HP equation, and the gammas from this LP, and over all pts

% INPUTS:
%  Ain: the data matrix
%  qparams: struct containing control parameters:
%    .q: MIO minimizes the qth deviation from the hyperplane. 
%      - if q > 0 it is the percentile of the m points in the data set
%      - if q < 0 then -q is the actual ordinal number (not percentile)
%      - if q = 0 then take the q from CBgen
%    .maxDist: points closer than this distance to a hyperplane are
%      "close". Distance is Euclidean.
%      There are 2 cases (0 is an error):
%         < 0: -maxDist is the percentile of distances from the first
%               hyperplane. Note it is in %. 
%               NOTE: maxDist = -16 IS HIGHLY RECOMMENDED.
%               If maxDist is not specified, then -16 is used.
%         > 0: an actual Euclidean distance to define maxDist. 
% OUTPUTS:
%  output: the results for CBqgen
%    .q: the integer value associated with the input qparams.q
%    .maxDist: the maximum Euclidean distance for "close" points as found
%              by CBgen
%    .CBgen time: time taken the run CBgen
%    .RHS: the solution hyperplane RHS constant
%    .weights: the solution hyperplane weights
%    .gammaLP: gamma for the solution hyperplane on the point subset
%    .gamma: gamma using solution hyperplane on all points.
%  inc: the output structure for CBgen
%
% DEPENDENCIES: this routine calls:
%   External solver Gurobi (free academic license).
%   CBgen (available on this Github site).
%   CBgen calls external solver MOSEK under certain conditions
%      (free academic license).

function [output,inc] = CBqgen(Ain,qparams)

tStart = tic;

% get m, n, and specific q
m = size(Ain,1);
n = size(Ain,2);

% set an integer value of q
q = qparams.q;
if q == 0
    fprintf("q = 0 is invalid. Aborting.\n")
    return
end
if q > 0
    % input q is a percentile
    q = floor(q/100*m);
else
    if q < 0
        % input q is a specific integer value
        q = -q;
    end
end
fprintf("q set at %d\n",q)
output.q = q;

%--------------------------------------------------------------------------
% STEP 1: run CBgen

fprintf("-----CBgen starts-----\n")
CBparam.maxDist = qparams.maxDist;
CBparam.mgood = 0;
[inc] = CBgen(Ain,CBparam);
fprintf("-----CBgen ends-----\n")
maxDist = inc.maxDist;
output.maxDist = maxDist;
fprintf("  CBgen sets maxDist at %f\n",maxDist)

output.CBgenTime = toc(tStart);

if inc.status == -1
    fprintf("  CBgen failure: aborting solution.\n")
    output.status = "ABORTED";
    output.gamma = Inf;
    output.gammaLP = Inf;
    return
end
if inc.status == 1
    fprintf("  CBgen exact solution: aborting solution: gamma is zero.\n")
    output.status = "Exact solution";
    output.gamma = 0.0;
    output.gammaLP = 0.0;
    return
end

%--------------------------------------------------------------------------
% STEP 2: select the q closest points based on the CBgen output hyperplane

% Calculate the point distances from the hyperplane
dist = Ain*inc.weightsOut - inc.RHSOut;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(inc.weightsOut);
edist = abs(dist/gradLen);
sortedEdist = [(1:m)',edist];
sortedEdist = sortrows(sortedEdist,2);
% Construct subset of data using only the q closest points
B = zeros(q,n);
% z indicates which rows are among the q closest
z = zeros(m,1);
for i=1:q
    z(sortedEdist(i,1),1) = 1;
    B(i,:) = Ain(sortedEdist(i,1),:);
end

%--------------------------------------------------------------------------
% STEP 3: solve LP to minimize maximum error on B

% Set up the gurobi parameters
gbparams = struct();
gbparams.OutputFlag = 0;

% Initialize the constraint matrix
% row order: q elastic constraints, q error constraints
% col order: n w, q eplus, q eminus, 1 emax
model.A = sparse(2*q,n+2*q+1);
model.rhs = zeros(2*q,1);
% q elastic cons, of form w1x1 + w2x2 + ... +wnxn + eplus - eminus = RHS. 
% Note use of the RHS output by CBgen
output.RHS = inc.RHSOut;
model.A(1:q,1:n+2*q) = [sparse(B),speye(q),-speye(q)];
model.rhs(1:q,1) = zeros(q,1) + inc.RHSOut;
model.sense(1:q,1) = repmat('=',q,1);
% q error constraints of form eplus_i + eminus_i - emax <= 0
model.A(q+1:2*q,n+1:n+2*q) = [speye(q),speye(q)];
model.A(q+1:2*q,n+2*q+1) = sparse(-ones(q,1));
model.sense(q+1:2*q,1) = repmat('<',q,1);
% variable lower bounds
model.lb = [zeros(n,1)-Inf;zeros(q*2+1,1)];
% variable upper bounds
model.ub = zeros(n+2*q+1,1)+Inf;
% specify variable types
model.vtype = repmat('C',n+2*q+1,1);
% set up objective function to minimize gamma
model.obj = [zeros(n+2*q,1);1.0];
model.modelsense = 'min';

% solve LP to minimize emax and return the HP equation
result = gurobi(model, gbparams);
output.weights = result.x(1:n,1);
output.RHS = inc.RHSOut;
output.gammaLP = result.x(n+2*q+1,1);

% calculate gamma over all points relative to the output hyperplane
% Calculate the point distances from the hyperplane
edist = abs(Ain*output.weights - output.RHS);
sortedEdist = sort(edist);
output.gamma = sortedEdist(q,1);

return
end

