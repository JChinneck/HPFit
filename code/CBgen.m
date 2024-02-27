% December 14, 2023
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA

% This function fits a hyperplane to a set of data points in a way that
% heuristically maximizes the number of points that are "close" to it. It
% uses a 4 step process: (1) fit a first hyperplane to the complete data
% set (and perhaps use it to caculate the maximum distance to a "close"
% point), (2) remove outliers using a new robust process and find a 2nd
% hyperplane, (3) reinstate any points that are not outliers relative to
% the 2nd hyperplane and remove any that are and place a 3rd hyperplane, 
% (4) repeat the last step for outliers relative to the 3rd hyperplane and 
% place the 4th and final hyperplane.
%    The Absolute Better Method may override the final hyperplane choice.
% INPUTS:
%   Aorig: the original data matrix (points x variables)
%   inParam: input parameters:
%     .mgood: if > 0 this is the known number of inlier points, which
%            must be the first mgood points in the set. Knowing this allows
%            the calculation of useful statistics for testing purposes. Of
%            course mgood is not known in general. Must be specified.
%     .maxDist: points closer than this distance to a hyperplane are
%              "close". Distance is Euclidean.
%              There are 3 cases:
%              < 0: -maxDist is the percentile of distances from the first
%                   hyperplane. Note it is in %. 
%                   NOTE: maxDist = -50 IS HIGHLY RECOMMENDED.
%                   If maxDist is not specified, then -50 is used.
%              = 0: means that closeAll is not used to help identify the
%                   best hyperplane. Result is just the 3rd hyperplane.
%              > 0: an actual Euclidean distance to define maxDist.
% OUTPUTS: these are all fields of inc. [x] has values 1,2,3,4,Out
%     .status: -1: failure, 0: OK, 1: exact solution (usually because m=n)
%     .maxDist: the final value of maxDist, points closer than maxDist to
%       a hyperplane are counted as being "close" to it 
%     .m[x]: the number of data points when finding the hyperplanes
%       [NOTE: hyperplane equations have this form:
%       weight_1*x_1 + weight_2*x_2 + .... = RHS]
%    .weights[x]: the weights in the hyperplane equations
%    .RHS[x]: the right hand side values in the hyperplane equations
%    .q: number of inliers estimated by the outFinder routine
%    .qout: number of inliers estimated by final hyperplane
%    .totSqDistAll[x]: total squared distance to all points
%    .solTime: solution time in seconds.
%    If maxDist ~= 0: 
%      .inc.closeAll[x]: number of points "close" to the hyperplane
%      If mgood > 0:
%        .closeTru[x]: number of mgood points "close" to the hyperplane
%    if mgood > 0:
%      .totSqDistTru[x]: total squared distances to the mgood points
%      .TSEstar[x]: the sum of mgood smallest squared errors for the 
%        hyperplanes. This is the same as "trimmed squared error" for 
%        mgood points. Note that it isn't necessarily the mgood non-outlier 
%        points identified on input, just mgood points in total.
%      .bnd2: this is an upper bound on the best possible TSEstar value,
%        calculated by fitting the PCA to only the mgood points and then
%        calculating the TSEstar for those points. 

function [inc] = CBgen2023(Aorig,inParam)
inc.status = 0;

fprintf("Input parameters:\n");
% mgood is the number of inliers, if known (useful for testing purposes)
if isfield(inParam,'mgood') == 1
    mgood = inParam.mgood;
else
    mgood = 0;
end
fprintf("  mgood %d\n",inParam.mgood)
if isfield(inParam,'maxDist') == 1
    maxDist = inParam.maxDist;
else
    maxDist = -50;
end
fprintf("  maxDist %f\n",inParam.maxDist)

tic;

% Get data table dimensions
m = size(Aorig,1);
norig = size(Aorig,2);

if m < norig
    % Too few points: abort
    fprintf("  Too few points: aborting. m = %d, n = %d.\n",m,norig)
    inc.weightsOut = zeros(norig,1);
    inc.RHSOut = 0;
    inc.edistOut = zeros(m,1);
    inc.bnd2 = 0;
    inc.status = -1;
    return
end

% Print some stats to console
if mgood > 0
    fprintf("  Stats: mgood %d mtot %d n %d mout %d mtot/n %f outFrac %f\n",...
        mgood,m,norig,m-mgood,m/norig,(m-mgood)/m)
else
    fprintf("  Stats: mtot %d n %d mtot/n %f\n",m,norig,m/norig)
end

%--------------------------------------------------------------------------
% initial PCA solution for the original dataset (HP1)

fprintf("Initial PCA on all %d points.\n",m)
inc.m1 = m;
[weights,RHS] = getPCAHP(Aorig);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

if max(edist) < 1.0e-6
    % Exact fit, within tolerance
    fprintf("  Exact fit. Max absolute Euclidean error: %f. Note m = %d and n = %d.\n",...
        max(edist),m,norig)
    inc.weightsOut = w;
    inc.RHSOut = w0;
    inc.edistOut = edist;
    inc.bnd2 = 0;
    inc.status = 1;
    inc.outStep = 1;
    return
end

% Find maxDist automatically if not prespecified
if maxDist < 0
    maxDist = prctile(edist,-maxDist);
    fprintf("  maxDist automatically selected as %f\n",maxDist)
end
inc.maxDist = maxDist; 

% Gather statistics about the initial PCA solution
inc.totSqDistAll1 = norm(edist(:,1).*edist(:,1),1);
inc.totSqDistAllOut = inc.totSqDistAll1;
if maxDist ~= 0
    inc.closeAll1 = sum(edist <= maxDist);
    inc.closeAllOut = inc.closeAll1;
    fprintf("  %d close points\n",inc.closeAll1)
end
inc.weights1 = weights;
inc.RHS1 = RHS;
inc.weightsOut = inc.weights1;
inc.RHSOut = inc.RHS1;

if mgood > 0
    inc.totSqDistTru1 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    inc.totSqDistTruOut = inc.totSqDistTru1;
    if maxDist ~= 0
        inc.closeTru1 = sum(edist(1:mgood,1) <= maxDist);
    end
    % Calculate TSEstar, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.TSEstar1 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    inc.TSEstarOut = inc.TSEstar1;
    fprintf("  TSEstar %f\n",inc.TSEstar1)
    %Calculate bnd2
    [weights,RHS] = getPCAHP(Aorig(1:mgood,:));
    dist = Aorig*weights - RHS;
    % Calculate the Euclidean point distances from the hyperplane
    gradLen = norm(weights);
    edist = abs(dist/gradLen);
    inc.bnd2 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    fprintf("  bnd2: %f\n",inc.bnd2)
end
inc.outStep = 1;

%--------------------------------------------------------------------------
% Analyze and remove outliers, then create HP2

% outFinder removes no more than mtot-n points, so that PCA can run
[OM] = outFinder(Aorig,mgood);
inc.q = OM.q;
inc.numZeroMedians = OM.numZeroMedians;
B = zeros(m,norig);
icount = 0;
for i=1:m
    if OM.TF(i,1) < 1
        icount = icount + 1;
        B(icount,:) = Aorig(i,:);
    end
end
B = B(1:icount,:);
inc.m2 = icount;
fprintf("Outliers removed: %d points remain.\n",icount)

% Get the PCA solution
[weights,RHS] = getPCAHP(B);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

inc.totSqDistAll2 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll2 = sum(edist <= maxDist);
    fprintf("  %d close points\n",inc.closeAll2)
end
inc.weights2 = weights;
inc.RHS2 = RHS;

if mgood > 0
    fprintf("  Fractions removed: tru %f out %f.\n",OM.outTruFrac,OM.outOutFrac)
    inc.outTruFrac = OM.outTruFrac;
    inc.outOutFrac = OM.outOutFrac;
    inc.totSqDistTru2 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru2 = sum(edist(1:mgood,1) <= maxDist);
    end
    % Calculate TSEstar, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.TSEstar2 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    fprintf("  TSEstar %f\n",inc.TSEstar2)
end

% Update output solution if appropriate
if maxDist ~= 0
    if inc.closeAll2 > inc.closeAll1
        inc.outStep = 2;
        inc.closeAllOut = inc.closeAll2;
        inc.weightsOut = inc.weights2;
        inc.RHSOut = inc.RHS2;
        inc.totSqDistAllOut = inc.totSqDistAll2;
        if mgood > 0
            inc.TSEstarOut = inc.TSEstar2;
            inc.totSqDistTruOut = inc.totSqDistTru2;
            inc.closeTruOut = inc.closeTru2;
        end
    end
end

%--------------------------------------------------------------------------
% Reinstate any nonoutliers and find HP3

% Too many points are sometimes removed for PCA to run, so guard against
outliers = isoutlier(edist);
if m - sum(outliers) < norig
    outliers = zeros(m,1);
    sortededist = sortrows([(1:m)',edist],2,'descend');
    for i = 1:(m-norig)
        outliers(sortededist(i,1),1) = 1;
    end
end

% Calculate outlier fractions at this stage
if mgood > 0
    inc.HP3FracInOut = sum(outliers(1:mgood)/mgood);
    inc.HP3FracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end
   
B = zeros(m,norig);
mB = 0;
for i=1:m
    if ~outliers(i,1)
        mB = mB + 1;
        B(mB,:) = Aorig(i,:);
    end
end
B = B(1:mB,:);
inc.m3 = mB;
fprintf("Reinstatement phase: %d points\n",mB)

% Get the PCA solution
[weights,RHS] = getPCAHP(B);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

inc.totSqDistAll3 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll3 = sum(edist <= maxDist);
    fprintf("  %d close points\n",inc.closeAll3)
end
inc.weights3 = weights;
inc.RHS3 = RHS;

if mgood > 0
    inc.totSqDistTru3 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru3 = sum(edist(1:mgood,1) <= maxDist);
    end
    % Calculate TSEstar, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.TSEstar3 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    fprintf("  TSEstar %f\n",inc.TSEstar3)
end

% Update output solution if appropriate
if maxDist ~= 0
    if inc.closeAll3 >= max(inc.closeAll1,inc.closeAll2)
        inc.outStep = 3;
        inc.closeAllOut = inc.closeAll3;
        inc.weightsOut = inc.weights3;
        inc.RHSOut = inc.RHS3;
        inc.totSqDistAllOut = inc.totSqDistAll3;
        if mgood > 0
            inc.TSEstarOut = inc.TSEstar3;
            inc.totSqDistTruOut = inc.totSqDistTru3;
            inc.closeTruOut = inc.closeTru3;
        end
    end
end

if maxDist == 0
    inc.outStep = 3;
    inc.weightsOut = inc.weights3;
    inc.RHSOut = inc.RHS3;
    inc.totSqDistAllOut = inc.totSqDistAll3;
    if mgood > 0
        inc.TSEstarOut = inc.TSEstar3;
        inc.totSqDistTruOut = inc.totSqDistTru3;
    end
end

%--------------------------------------------------------------------------
% Check phase. Remove outliers again and find HP4.

% Avoid removing too many points. We need m for PCA to run
outliers = isoutlier(edist);
if m - sum(outliers) < norig
    outliers = zeros(m,1);
    sortededist = sortrows([(1:m)',edist],2,'descend');
    for i = 1:(m-norig)
        outliers(sortededist(i,1),1) = 1;
    end
end

% Calculate outlier fractions at this stage
if mgood > 0
    inc.HP4FracInOut = sum(outliers(1:mgood)/mgood);
    inc.HP4FracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end
   
B = zeros(m,norig);
mB = 0;
for i=1:m
    if ~outliers(i,1)
        mB = mB + 1;
        B(mB,:) = Aorig(i,:);
    end
end
B = B(1:mB,:);
inc.m4 = mB;
fprintf("Check phase: %d points\n",mB)

% Get the PCA solution
[weights,RHS] = getPCAHP(B);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

inc.totSqDistAll4 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll4 = sum(edist <= maxDist);
    fprintf("  %d close points\n",inc.closeAll4)
end
inc.weights4 = weights;
inc.RHS4 = RHS;

if mgood > 0
    inc.totSqDistTru4 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru4 = sum(edist(1:mgood,1) <= maxDist);
    end
    % Calculate TSEstar, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.TSEstar4 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    fprintf("  TSEstar %f\n",inc.TSEstar4)
end

% Update output solution if appropriate
if maxDist ~= 0
    if inc.closeAll4 >= inc.closeAllOut
        inc.outStep = 4;
        inc.closeAllOut = inc.closeAll4;
        inc.weightsOut = inc.weights4;
        inc.RHSOut = inc.RHS4;
        inc.totSqDistAllOut = inc.totSqDistAll4;
        if mgood > 0
            inc.TSEstarOut = inc.TSEstar4;
            inc.totSqDistTruOut = inc.totSqDistTru4;
            inc.closeTruOut = inc.closeTru4;
        end
    end
end

if maxDist == 0
    inc.outStep = 4;
    inc.weightsOut = inc.weights4;
    inc.RHSOut = inc.RHS4;
    inc.totSqDistAllOut = inc.totSqDistAll4;
    if mgood > 0
        inc.TSEstarOut = inc.TSEstar4;
        inc.totSqDistTruOut = inc.totSqDistTru4;
    end
end

%--------------------------------------------------------------------------
% Collate and print final solution 

inc.solTime = toc;

fprintf("Final solution from HP %d.",inc.outStep)
if inc.outStep ~= 4
    fprintf(" ABM active.\n")
else
    fprintf(" ABM inactive.\n")
end
if maxDist ~= 0
    fprintf("  %d close points\n",inc.closeAllOut)
end
if mgood > 0
    fprintf("  TSEstar %f\n",inc.TSEstarOut)
end

% Calculate the estimated number of inliers at output (inc.qout)
% Calculate the point distances from the hyperplane
dist = Aorig*inc.weightsOut - inc.RHSOut;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(inc.weightsOut);
edist = abs(dist/gradLen);
outliers = isoutlier(edist);
inc.qout = m - sum(outliers);
% Calculate outlier fractions at final fit
if mgood > 0
    inc.FinalFracInOut = sum(outliers(1:mgood)/mgood);
    inc.FinalFracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end

return
end
% =========================================================================
% Analyzes the input data set to identify outliers.
% Looks along the columns as well as along the columns of the principal
% components. Finds the largest relative distance from the median distance
% to the median for each point in the 2n columns mentioned. Uses an offset
% in case the median (used as a denominator) is zero, or close to it.
% if mgood > 0 then it returns measures on the accuracy of identifying
% outliers.
% INPUTS:
%   A: data matrix
%   mgood: number of non-outliers, if known, in which case you get
%          statistics are generated on how accurate the outlier
%          identification is. If mgood <= 0, or mgood >= m, then no such
%          statistics are generated. This is for testing purposes.
% OUTPUTS: these are all fields in outMeasure
%   .max: the maximum relative distance from the median, in medians, in any
%         axis, including the original variables and the PCA axes
%   .TF: mx1 logical vector listing points identified as outliers
%   .count: total number of points identified as outliers
%   .q: the estimated number of inliers
%   If mgood > 0 and mgood < m:
%     .outTru: number of nonoutlier points identified as outliers
%     .outTruFrac: fraction of nonoutlier points identified as outliers
%     .outOut: number of outlier points identified as outliers
%     .outOutFrac: fraction of outlier points identified as outliers
function [outMeasure] = outFinder(A,mgood)
m = size(A,1);
n = size(A,2);
outMeasure.numZeroMedians = 0;
outMeasure.max = zeros(m,1) - Inf;

% Axis-aligned hyperplanes
for j=1:n
    absDiffs = abs(A(:,j) - median(A(:,j)));
    if median(absDiffs) < 1.0e-6
        outMeasure.numZeroMedians = outMeasure.numZeroMedians + 1;
    end
    % Offset in case of zero median
    fracDiffs = absDiffs/(median(absDiffs)+1);
    for i=1:m
        if fracDiffs(i,1) > outMeasure.max(i,1)
            outMeasure.max(i,1) = fracDiffs(i,1);
        end
    end
end

% principal component axes
[~,score,~] = pca(A);
for j=1:n
    absDiffs = abs(score(:,j) - median(score(:,j)));
    if median(absDiffs) < 1.0e-6
        outMeasure.numZeroMedians = outMeasure.numZeroMedians + 1;
    end
    % Offset in case of zero median
    fracDiffs = absDiffs/(median(absDiffs)+1);
    for i=1:m
        if fracDiffs(i,1) > outMeasure.max(i,1)
            outMeasure.max(i,1) = fracDiffs(i,1);
        end
    end
end

% Use the first abrupt change value in sorted outMeasure.max after the 
% starting point to identify outliers
sortedMax = sort(outMeasure.max);
changes = ischange(sortedMax,'linear');
% Set initial defaults for the cutoff value. If there is no abrupt change
% then the default cutoff values are used.
if m/n <= 2
    % Can't start at m/2 or we don't leave enough points for PCA to run
    istart = n + 1;
    cutoff = sortedMax(istart,1);
else
    istart = ceil(m/2);
    cutoff = sortedMax(istart,1);
end
outMeasure.q = istart - 1;
 
for i = istart:m
    if changes(i,1)
        cutoff = sortedMax(i,1);
        outMeasure.q = i-1;
        break
    end
end
outMeasure.TF = outMeasure.max >= cutoff;
outMeasure.count = sum(outMeasure.TF);

% if mgood is nonzero, then calculate some accuracy measures
if (mgood > 0) && (mgood < m) 
    outMeasure.outTru = sum(outMeasure.TF(1:mgood,1));
    outMeasure.outTruFrac = outMeasure.outTru/mgood;
    outMeasure.outOut = sum(outMeasure.TF(mgood+1:m,1));
    outMeasure.outOutFrac = outMeasure.outOut/(m-mgood);
end

return
end
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
