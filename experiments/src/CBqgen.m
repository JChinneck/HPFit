% August 4, 2022
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% Heuristic function to minimize the qth error in a general data set. Main
% steps are:
%   1. Run CBgen to obtain hyperplane fit.
%   2. Select the q points closest to the the CBgen HP
%   3. Cycle until no reduction in gamma:
%      3.1 Solve LP on the q retained points to minimize the maximum error
%      3.2 Calculate gamma from resulting HP
%      3.3 Select the q points closest to HP

% INPUTS:
%  Ain: the data matrix
%  qparams: struct containing control parameters:
%    .q: MIO minimizes the qth deviation from the hyperplane. 
%      - if q > 0 it is the percentile of the m points in the data set
%      - if q < 0 then -q is the actual ordinal number (not percentile)
%      - q = 0 is an error
%    .maxDist: points closer than this distance to a hyperplane are
%      "close". Distance is Euclidean.
%      There are 3 cases:
%         < 0: -maxDist is the percentile of distances from the first
%               hyperplane. Note it is in %. 
%               NOTE: maxDist = -16 IS HIGHLY RECOMMENDED.
%               If maxDist is not specified, then -16 is used.
%         = 0: means that closeAll is not used to help identify the
%              best hyperplane. Result is just the final hyperplane.
%         > 0: a Euclidean distance value for maxDist. 
% OUTPUTS:
%  output: the results for CBqgen
%    .status: CBgen exit status
%    .q: the integer value associated with the input qparams.q
%    .maxDist: the maximum Euclidean distance for "close" points as found
%              by CBgen
%    .CBgen time: time taken the run CBgen
%    .RHS: the solution hyperplane RHS constant
%    .weights: the solution hyperplane weights
%    .gammaLP: gamma for the solution hyperplane on the point subset
%    .gamma: gamma using solution hyperplane on all points.
%    .gammaN: this is .gamma scaled to a right-hand-side constant of n, the
%       number of variables, for a fair comparison with other methods
%    .z: binary vector indicating the q closest points
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
    output.status = "Aborted_q=0";
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
output.status = "CBgen_succeeds";
if inc.status == -1
    fprintf("  CBgen failure: aborting solution.\n")
    output.status = "CBgen_aborted";
    output.gamma = Inf;
    output.gammaLP = Inf;
    return
end
if inc.status == 1
    fprintf("  CBgen exact solution: aborting MIO solution.\n")
    output.status = "CBgen_exact_solution";
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
for i=1:q
    B(i,:) = Ain(sortedEdist(i,1),:);
end

%--------------------------------------------------------------------------
% STEP 3: solve LP to minimize maximum error on B, in a cycle

% Set up the gurobi parameters
gbparams = struct();
gbparams.OutputFlag = 0;
gbparams.Threads=14;

gammaBest = Inf;
for itn = 1:m
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
    weights = result.x(1:n,1);
    RHS = inc.RHSOut;
    gammaLP = result.x(n+2*q+1,1);
    
    % calculate gamma over all points relative to the output hyperplane
    % Calculate the absolute point distances from the hyperplane
    dist = abs(Ain*weights - RHS);
    sortedDist = [(1:m)',dist];
    sortedDist = sortrows(sortedDist,2);
    gamma = sortedDist(q,2);
    gammaN = abs(gamma/RHS*n);
    
    if gammaN >= gammaBest
        break
    end
    gammaBest = gammaN;
    
    output.gamma = gamma;
    output.gammaN = gammaN;
    output.gammaLP = gammaLP;
    output.weights = weights;
    output.RHS = RHS;
    fprintf ("Gamma is %f at itn %d.\n",output.gammaN,itn)

    % Set up for the next iteration
    % output.z indicates the q points having the smallest error
    B = zeros(q,n);
    output.z = zeros(m,1);
    for i=1:q
        output.z(sortedDist(i,1),1) = 1;
        B(i,:) = Ain(sortedDist(i,1),:);
    end
end
return
end

