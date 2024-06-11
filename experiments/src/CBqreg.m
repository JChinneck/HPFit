% August 4, 2022
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% Heuristic function to minimize the qth error in a regression data set. 
% Main steps are:
%   1. Run CBreg to obtain hyperplane fit.
%   2. Select the q points closest to the the CBreg HP
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
%    .maxResid: points closer than this distance to a hyperplane are
%      "close". Distance is Euclidean.
%      There are 2 cases (0 is an error):
%         < 0: -maxResid is the percentile of distances from the first
%               hyperplane. Note it is in %. 
%               NOTE: maxResid = -16 IS HIGHLY RECOMMENDED.
%               If maxDist is not specified, then -16 is used.
%         = 0: means that closeAll is not used to help identify the
%              best hyperplane. Result is just the final hyperplane.
%         > 0: a value for maxResid. 
% OUTPUTS:
%  output: the results for CBqreg
%    .status: CBreg exit status
%    .q: the integer value associated with the input qparams.q
%    .maxResid: the maximum error for "close" points as found by CBreg
%    .CBreg time: time taken to run CBreg
%    .RHS: the solution hyperplane RHS constant
%    .w: the solution hyperplane weights
%    .w0: the solution hypeplane constant
%    .gammaLP: gamma for the solution hyperplane on the point subset
%    .gamma: gamma using solution hyperplane on all points.
%    .z: binary vector indicating the q points with the smallest residuals
%  inc: the output structure for CBreg
%
% DEPENDENCIES: this routine calls:
%   External solver Gurobi (free academic license).
%   CBgen (available on this Github site).
%   CBgen calls external solver MOSEK under certain conditions
%      (free academic license).

function [output,inc] = CBqreg(y,Ain,qparams)

tStart = tic;

% The regression equation is w0 + w1x1 + w2x2 + ... + wnxn = y,
% so a column of ones is added at the beginning of Ain for w0
A = [ones(size(Ain,1),1),Ain];
m = size(A,1);
n = size(A,2);

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
% STEP 1: run CBreg

fprintf("-----CBreg starts-----\n")
CBparam.maxResid = qparams.maxResid;
CBparam.mgood = 0;
[inc] = CBreg(y,Ain,CBparam);
fprintf("-----CBreg ends-----\n")
maxResid = inc.maxResid;
output.maxResid = maxResid;
fprintf("  CBreg sets maxResid at %f\n",maxResid)

output.CBregTime = toc(tStart)
output.status = "CBreg_succeeds";
if inc.status == -1
    fprintf("  CBreg failure: aborting solution.\n")
    output.status = "CBreg_aborted";
    output.gamma = Inf;
    output.gammaLP = Inf;
    return
end
if inc.status == 1
    fprintf("  CBreg exact solution: aborting solution. Gamma = 0.\n")
    output.status = "CBreg_exact_solution";
    output.gamma = 0.0;
    output.gammaLP = 0.0;
    return
end

%--------------------------------------------------------------------------
% STEP 2: select the q closest points based on the CBreg output hyperplane

% Calculate the errors relative to the CBreg hyperplane
absResid = abs(inc.residOut);
sortedResid = [(1:m)',absResid];
sortedResid = sortrows(sortedResid,2);
% Construct subset of data using only the q closest points
B = zeros(q,n);
yB = zeros(q,1);
for i=1:q
    B(i,:) = A(sortedResid(i,1),:);
    yB(i,1) = y(sortedResid(i,1),:);
end

%--------------------------------------------------------------------------
% STEP 3: solve LP to minimize maximum error on B, in a cycle.

% Set up the gurobi parameters
gbparams = struct();
gbparams.OutputFlag = 0; 
gbparams.Threads=14;
gammaBest = Inf;

for itn =1:m
    fprintf("iteration %d\n", itn);
    toc(tStart)
    % Initialize the constraint matrix
    % row order: q elastic constraints, q error constraints
    % col order: n w, q eplus, q eminus, 1 emax
    model.A = sparse(2*q,n+2*q+1);
    model.rhs = zeros(2*q,1);
    % q elastic cons, of form wo + w1x1 + w2x2 + ... +wnxn + eplus - eminus = y.
    model.A(1:q,1:n+2*q) = [sparse(B),speye(q),-speye(q)];
    model.rhs(1:q,1) = yB;
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
    
    w0 = result.x(1,1);
    w = result.x(2:n,1);
    gammaLP = result.x(n+2*q+1,1);
    
    % calculate gamma over all points relative to the output hyperplane
    % Calculate the point distances from the hyperplane
    absResid = abs(w0 + Ain*w - y);
    sortedAbsResid = [(1:m)',absResid];
    sortedAbsResid = sortrows(sortedAbsResid,2);
    gamma = sortedAbsResid(q,2);
    
    if gamma >= gammaBest
        break
    end
    gammaBest = gamma;
    
    output.gamma = gamma;
    output.w0 = w0;
    output.w = w;
    output.gammaLP = gammaLP;
    
    fprintf ("Gamma is %f at itn %d.\n",output.gamma,itn)
    
    % Set up for the next iteration
    % z marks the q points having smallest residuals
    output.z = zeros(m,1);
    B = zeros(q,n);
    yB = zeros(q,1);
    for i=1:q
        output.z(sortedAbsResid(i,1),1) = 1;
        B(i,:) = A(sortedAbsResid(i,1),:);
        yB(i,1) = y(sortedAbsResid(i,1),1);
    end
end
return
end
