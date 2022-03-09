% March 8, 2022
% John W. Chinneck, Systems and Computer Engineering, Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, Virginia Commonwealth University, Richmond, Virginia, USA

% This function fits a hyperplane to a set of data points in a way that
% heuristically maximizes the number of points that are "close" to it. It
% uses a 3 step process: (1) fit a first hyperplane to the complete data
% set (and perhaps use it to caculate the maximum distance to a "close"
% point), (2) remove outliers using a new process and find a 2nd
% hyperplane, (3) reinstate any points that are not outliers relative to
% the 2nd hyperplane and place a 3rd and final hyperplane.
% INPUTS:
%   Aorig: the original data matrix (points x variables)
%   inParam: input parameters:
%     mgood: if > 0 this is the known number of nonoutlier points, which
%            must be the first mgood points in the set. Knowing this allows
%            the calculation of useful statistics for testing purposes. Of
%            course mgood is not known in general.
%     maxDist: points closer than this distance to a hyperplane are
%              "close". Distance is Euclidean.
%              There are 3 cases:
%              < 0: -maxDist is the percentile of distances from the first
%                   hyperplane to used as maxDist. Note it is in %. 
%                   maxDist = -16 is recommended.
%              = 0: means that maxDist is not used to help identify the
%                   best hyperplane to output.
%              > 0: used as an actual Euclidean distance to define maxDist
%                   which defines close points.
% OUTPUTS: these are all fields of inc
%     .maxDist: the final value of maxDist, points closer than maxDist to
%       a hyperplane are counted as being "close" to it 
%     .m1, .m2, .m3: the number of data points when finding the three
%       hyperplanes
%    [NOTE: hyerplane equations have this form:
%       weight_1*x_1 + weight_2*x_2 + .... = RHS]
%    .weights1, .weights2, .weights3, .weightsOut: the weights in the three
%       hyperplane equations, and the weights in the output hyperplane
%    .RHS1, .RHS2, .RHS3, .RHSOut: the right hand side values in the three
%       hyperplane equations, and the right hand side in the output
%       hyperplane equation.
%    .totSqDistAll1, totSqDistAll2, totSqDistAll3, totSqDistAllOut: total
%       squared distance to all points
%    .solTime: solution time in seconds.
%    If maxDist ~= 0: 
%      .inc.closeAll1, inc.closeAll2, inc.closeAll3, inc.closeAllOut: 
%        number of points "close" to the hyperplane
%      If mgood > 0:
%        .closeTru1, .closeTru2, .closeTru3, .closeTruOut:
%          the number of points "close" to the hyperplane
%    if mgood > 0:
%      .totSqDistTru1, .totSqDistTru2, .totSqDistTru3, .totSqDistTruOut:
%        the total squared distances to the mgood points
%      .SMSSE1, .SMSSE2, .SMSSE3, .SMSSEOut: the "sum of mgood smallest
%        squared errors" for the hyperplanes. This is the same as "trimmed 
%        squared error" for mgood points. Note that isn't necessarily the
%        mgood non-outlier points identified as input, just mgood points in
%        total.
%      .bnd2: this is a tight upper bound on the best possible SMSSE value,
%        calculated by fitting the PCA to only the mgood points and then
%        calculating the SMSSE for those points. 

function [inc] = CBgen(Aorig,inParam)

mgood = inParam.mgood;
tic;

% Get data table dimensions
m = size(Aorig,1);
norig = size(Aorig,2);
% Print some stats to console
if mgood > 0
    fprintf("mgood %d mtot %d n %d mout %d mtot/n %f outFrac %f\n",mgood,m,norig,m-mgood,m/norig,(m-mgood)/m)
else
    fprintf("mtot %d n %d mtot/n %f\n",m,norig,m/norig)
end

% initial PCA solution for the original dataset ---------------------------
inc.m1 = m;
[weights,RHS] = getPCAHP(Aorig);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

% Find maxDist automatically if not prespecified
maxDist = inParam.maxDist;
if inParam.maxDist < 0
    maxDist = prctile(edist,-maxDist);
    fprintf("maxDist automatically selected as %f\n",maxDist)
end
inc.maxDist = maxDist; 

% Gather statistics about the initial PCA solution
inc.totSqDistAll1 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll1 = sum(edist <= maxDist);
end
inc.weights1 = weights;
inc.RHS1 = RHS;
if mgood > 0
    inc.totSqDistTru1 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru1 = sum(edist(1:mgood,1) < maxDist);
    end
    fprintf("Initial PCA: tot sq dist tru %f\n",inc.totSqDistTru1)  
    % Calculate SMSSE, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.SMSSE1 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    inc.SMSSEout = inc.SMSSE1;
    fprintf("SMSSE1 %f\n",inc.SMSSE1)
    %Calculate bnd2
    [weights,RHS] = getPCAHP(Aorig(1:mgood,:));
    dist = Aorig*weights - RHS;
    % Calculate the Euclidean point distances from the hyperplane
    gradLen = norm(weights);
    edist = abs(dist/gradLen);
    inc.bnd2 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    inc.SMSSEout = inc.SMSSE1;
end

% Analyze and remove outliers ---------------------------------------------
% outFinder removes no more than mtot-n points, so that PCA can run
[OM] = outFinder(Aorig,mgood);
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

% get the PCA solution
[weights,RHS] = getPCAHP(B);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

inc.totSqDistAll2 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll2 = sum(edist <= maxDist);
end
inc.weights2 = weights;
inc.RHS2 = RHS;

if mgood > 0
    fprintf("Fractions removed: tru %f out %f. %d total pts left.\n",OM.outTruFrac,OM.outOutFrac,icount)
    inc.outTruFrac = OM.outTruFrac;
    inc.outOutFrac = OM.outOutFrac;
    inc.totSqDistTru2 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru2 = sum(edist(1:mgood,1) < maxDist);
    end
    fprintf("After removal: tot sq dist tru %f\n",inc.totSqDistTru2)
    % Calculate SMSSE, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.SMSSE2 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    fprintf("SMSSE2 %f\n",inc.SMSSE2)
end

% Update output solution if appropriate
if maxDist ~= 0
    if inc.closeAll2 > inc.closeAll1
        inc.weightsOut = inc.weights2;
        inc.RHSOut = inc.RHS2;
        inc.totSqDistAllOut = inc.totSqDistAll2;
        if mgood > 0
            inc.SMSSEout = inc.SMSSE2;
            inc.totSqDistTruOut = inc.totSqDistTru2;
            inc.closeTruOut = inc.closeTru2;
        end
    end
end

% reinstate any nonoutliers and find HP again -----------------------------
% Too many points are sometimes removed for PCA to run, so guard against
outliers = isoutlier(edist);
if m - sum(outliers) < norig
    outliers = zeros(m,1);
    sortededist = sortrows([(1:m)',edist],2,'descend');
    for i = 1:(m-norig)
        outliers(sortededist(i,1),1) = 1;
    end
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
fprintf("%d points after reinstatement.\n",mB)

% get the PCA solution
[weights,RHS] = getPCAHP(B);
% Calculate the point distances from the hyperplane
dist = Aorig*weights - RHS;
% Calculate the Euclidean point distances from the hyperplane
gradLen = norm(weights);
edist = abs(dist/gradLen);

inc.totSqDistAll3 = norm(edist(:,1).*edist(:,1),1);
if maxDist ~= 0
    inc.closeAll3 = sum(edist <= maxDist);
end
inc.weights3 = weights;
inc.RHS3 = RHS;

if mgood > 0
    inc.totSqDistTru3 = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    if maxDist ~= 0
        inc.closeTru3 = sum(edist(1:mgood,1) < maxDist);
    end
    fprintf("After reinstatement: tot sq dist tru %f\n",inc.totSqDistTru3)
    % Calculate SMSSE, the sum of the mgood smallest squared errors
    sortededist = sort(edist);
    inc.SMSSE3 = norm(sortededist(1:mgood,1).*sortededist(1:mgood,1),1);
    fprintf("SMSSE3 %f\n",inc.SMSSE3)
end

% Update output solution if appropriate
if maxDist ~= 0
    if inc.closeAll3 >= max(inc.closeAll1,inc.closeAll2)
        inc.weightsOut = inc.weights3;
        inc.RHSOut = inc.RHS3;
        inc.totSqDistAllOut = inc.totSqDistAll3;
        if mgood > 0
            inc.SMSSEout = inc.SMSSE3;
            inc.totSqDistTruOut = inc.totSqDistTru3;
            inc.closeTruOut = inc.closeTru3;
        end
    end
end

if maxDist == 0
    inc.weightsOut = inc.weights3;
    inc.RHSout = inc.RHS3;
    inc.totSqDistAllOut = inc.totSqDistAll3;
    if mgood > 0
        inc.SMSSEout = inc.SMSSE3;
        inc.totSqDistTruOut = inc.totSqDistTru3;
    end
end

if mgood > 0
    fprintf("Final SMSSE %f\n",inc.SMSSEout)
end

inc.solTime = toc;

return
end
% -------------------------------------------------------------------------
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
%   If mgood > 0 and mgood < m:
%     .outTru: number of nonoutlier points identified as outliers
%     .outTruFrac: fraction of nonoutlier points identified as outliers
%     .outOut: number of outlier points identified as outliers
%     .outOutFrac: fraction of outlier points identified as outliers
function [outMeasure] = outFinder(A,mgood)
m = size(A,1);
n = size(A,2);
outMeasure.max = zeros(m,1) - Inf;

% Axis-aligned hyperplanes
for j=1:n
    absDiffs = abs(A(:,j) - median(A(:,j)));
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
if m/n < 2
    % Can't start at m/2 or we don't leave enough points for PCA to run
    istart = n + 1;
    cutoff = sortedMax(istart,1);
else
    istart = ceil(m/2);
    cutoff = sortedMax(istart,1);
end
 
for i = istart:m
    if changes(i,1)
        cutoff = sortedMax(i,1);
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
