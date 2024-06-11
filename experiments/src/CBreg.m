% December 14, 2023
% John W. Chinneck, Systems and Computer Engineering, 
%   Carleton University, Ottawa, Canada
% J. Paul Brooks, Dept. of Information Systems, 
%   Virginia Commonwealth University, Richmond, Virginia, USA
 
% This function fits a regression hyperplane to a set of data points to
% heuristically maximize the number of calculated points that are close to
% the given output points.
%   It uses a 4 step process: (1) regress on the complete data
% set (and perhaps use it to caculate the maximum residual to a "close"
% point), (2) remove outliers using a robust process and regress a 2nd
% hyperplane, (3) reinstate any points whose output residuals are not 
% outliers relative to the 2nd hyperplane, and remove any that are, and 
% regress a 3rd hyperplane, (4) repeat the last step relative to the 3rd 
% hyperplane and regress a 4th hyperplane. The Absolute Better Method may
% choose which of the 4 hyperplanes to output.
% INPUTS:
%   y: the output variable values
%   Aorig: the original data matrix (points x predictor variables)
%   inParam: input parameters:
%     .mgood: if > 0 this is the known number of nonoutlier points, which
%            must be the first mgood points in the y and Aorig. Knowing
%            them allows the calculation of useful statistics for testing 
%            purposes. Of course mgood is not known in general.
%            Must be specified on input.
%     .maxResid: points having a residual error smaller than this are
%              "close". 
%              There are 3 cases:
%              < 0: -maxResid is the percentile of residuals from the first
%                   regression to be used as maxDist. Note it is in %. 
%                   NOTE: maxDist = -50 IS HIGHLY RECOMMENDED.
%              = 0: means that maxResid is not used to help identify the
%                   best hyperplane to output.
%              > 0: used as an actual residual value to define maxResid
%                   which defines close points.
% OUTPUTS: these are all fields of inc. [x] has values 1,2,3,4,Out
%     .status: -1: failure, 0: OK, 1: exact solution (usually because m=n)
%     .maxResid: the final value of maxResid, points with residuals smaller
%       than this are counted as "close"
%     .m[x]: the number of data points when finding the three
%       regression hyperplanes
%    [NOTE: the regression hyerplane equations have this form:
%       w0 + w1x1 + w2x2 + ... + wnxn = y,]
%    .w0[x]: the constants in the hyperplane equations
%    .weights[x]: the weights in the hyperplane equations
%    .q: number of inliers estimated by the outFinder routine
%    .qout: number of inliers estimated by final hyperplane
%    .totSqResidAll[x]: total squared residual over all points
%    .solTime: solution time in seconds.
%    If maxDist ~= 0: 
%      .inc.numCloseAll[x]: number of "close" residuals
%      If mgood > 0:
%        .numCloseTru[x]: number of good points that have "close" residuals
%    If mgood > 0:
%      .totSqResidTru[x]: total squared residuals over the mgood points
%      .TSEstar[x]: the "sum of mgood smallest squared residual errors" for 
%        the hyperplanes. Note that isn't necessarily over the mgood 
%        pre-specified inlier points, just over mgood points in total.
%      .bnd2: this is an upper bound on the best possible TSEstar value,
%        calculated by fitting to only the mgood points and then
%        calculating the TSEstar for those points. Usually tight, but not
%        necessarily so high outlier fractions.
%      .RMSETru[x]: root mean square residual for the mgood points
%      .MSETru[x]: mean squared residual for the mgood points
%      .MAETru[x]: mean absolute residual for the mgood points

function [inc] = CBreg(y,Aorig,inParam)
inc.status = 0;

if isfield(inParam,'mgood') == 1
    mgood = inParam.mgood;
else
    mgood = 0;
end
if isfield(inParam,'maxResid') == 1
    maxResid = inParam.maxResid;
else
    % default to the recommended value
    maxResid = -50;
end
fprintf("Input parameters:\n")
fprintf("  maxResid %f\n",maxResid)
fprintf("  mgood %d\n",mgood)

tic;

% Get data table dimensions
m = size(Aorig,1);
norig = size(Aorig,2);

if m < norig
    % Too few points: abort
    fprintf("  Too few points: aborting. m = %d, n = %d.\n",m,norig)
    inc.weightsOut = zeros(norig,1);
    inc.w0Out = 0;
    inc.residOut = zeros(m,1);
    inc.maxResid = 0;
    inc.bnd2 = 0;
    inc.status = -1;
    return
end

if mgood > 0
    fprintf("  Stats: mgood %d mtot %d n %d mout %d mtot/n %f outFrac %f\n",mgood,m,norig,m-mgood,m/norig,(m-mgood)/m)
else
    fprintf("  Stats: mtot %d n %d mtot/n %f\n",m,norig,m/norig)
end

% The regression equation is w0 + w1x1 + w2x2 + ... + wnxn = y,
% so a column of ones is added at the beginning of Aorig
A = [ones(m,1),Aorig];

%--------------------------------------------------------------------------
% Initial regression solution for the original dataset

fprintf("Initial regression on all %d points.\n",m)
inc.m1 = m;
beta = regress(y,A);
w0 = beta(1,1);
w = beta(2:norig+1,1);
resid = w0 + Aorig*w - y;
absresid = abs(resid);
curr_time = toc;
fprintf("Initial regression complete, time %f\n", curr_time);

if max(absresid) < 1.0e-6
    % Exact fit, within tolerance
    fprintf("  Exact fit. Max absolute residual: %f. Note m = %d and n = %d.\n",...
        max(absresid),m,norig)
    inc.weightsOut = w;
    inc.w0Out = w0;
    inc.residOut = resid;
    inc.maxResid = 0;
    inc.bnd2 = 0;
    inc.status = 1;
    if mgood > 0
        inc.numCloseTruOut = mgood;
        inc.numCloseAllOut = mgood;
    end
    return
end

if maxResid < 0
    maxResid = prctile(absresid,-maxResid);
    fprintf("  maxResid automatically selected as %f\n",maxResid)
end
inc.maxResid = maxResid; % record the final value of maxResid

inc.resid1 = resid;
if mgood > 0
    if maxResid >= 0
        inc.numCloseTru1 = sum(absresid(1:mgood,1) <= maxResid);
        inc.numCloseTruOut = inc.numCloseTru1;
    end
    inc.totSqResidTru1 = norm(absresid(1:mgood,1).*absresid(1:mgood,1),1);
    inc.totSqResidTruOut = inc.totSqResidTru1;
    inc.RMSETru1 = sqrt(inc.totSqResidTru1);
    inc.MSETru1 = inc.totSqResidTru1/mgood;
    inc.MAETru1 = sum(abs(absresid(1:mgood,1)))/mgood;
    % Calculate bnd2: regression on only the mgood points
    betabnd = regress(y(1:mgood,1),A(1:mgood,:));
    w0bnd = betabnd(1,1);
    wbnd = betabnd(2:norig+1,1);
    residbnd = w0bnd + Aorig*wbnd - y;
    inc.bnd2 = norm(residbnd(1:mgood,1).*residbnd(1:mgood,1),1);
    fprintf("  bnd2: %f\n",inc.bnd2)
    % Calculate TSEstar
    sortedabsResid = sort(absresid);
    inc.TSEstar1 = norm(sortedabsResid(1:mgood,1).*sortedabsResid(1:mgood,1),1);
    inc.TSEstarOut = inc.TSEstar1;
    fprintf("  TSEstar %f. \n",inc.TSEstar1)
end
inc.totSqResidAll1 = norm(absresid(:,1).*absresid(:,1),1);
inc.weights1 = w;
inc.w01 = w0;
if maxResid > 0
    inc.numCloseAll1 = sum(absresid(:,1) <= maxResid);
    inc.numCloseAllOut = inc.numCloseAll1;
    fprintf("  %d close points\n",inc.numCloseAll1)
end

% Initialize the output values
outStep = 1;
inc.weightsOut = inc.weights1;
inc.w0Out = inc.w01;
inc.residOut = inc.resid1;
inc.totSqResidAllOut = inc.totSqResidAll1;

%--------------------------------------------------------------------------
% Analyze and remove outliers and place HP2

curr_time = toc;
fprintf("Starting outfinder %f\n", curr_time);
[OM] = outFinderReg(Aorig,y,mgood);

curr_time = toc;
fprintf("outfinder complete %f\n", curr_time);

inc.q = OM.q;
inc.numZeroMedians = OM.numZeroMedians;
B = zeros(m,norig);
y1 = zeros(m,1);
icount = 0;
for i=1:m
    if OM.TF(i,1) < 1
        icount = icount + 1;
        B(icount,:) = Aorig(i,:);
        y1(icount,1) = y(i,1);
    end
end
B = B(1:icount,:);
y1 = y1(1:icount,1);
inc.m2 = icount;

fprintf("Outliers removed: %d pts remain.\n",icount)
if mgood > 0
    fprintf("  Fractions removed: tru %f out %f\n",OM.outTruFrac,OM.outOutFrac)
    inc.outTruFrac = OM.outTruFrac;
    inc.outOutFrac = OM.outOutFrac;
end
curr_time = toc;
fprintf("removed outliers %f\n", curr_time);

% Get the regression solution
A = [ones(icount,1),B];
beta = regress(y1,A);
w0 = beta(1,1);
w = beta(2:norig+1,1);
resid = w0 + Aorig*w - y;
absresid = abs(resid);

curr_time = toc;
fprintf("second regression complete %f\n", curr_time);

inc.resid2 = resid;
if mgood > 0
    inc.totSqResidTru2 = norm(absresid(1:mgood,1).*absresid(1:mgood,1),1);
    inc.RMSETru2 = sqrt(inc.totSqResidTru2);
    inc.MSETru2 = inc.totSqResidTru2/mgood;
    inc.MAETru2 = sum(abs(absresid(1:mgood,1)))/mgood;
    if maxResid > 0
        inc.numCloseTru2 = sum(absresid(1:mgood,1) <= maxResid);
    end
    % Calculate TSEstar
    sortedabsResid = sort(absresid);
    inc.TSEstar2 = norm(sortedabsResid(1:mgood,1).*sortedabsResid(1:mgood,1),1);
    fprintf("  TSEstar %f.\n",inc.TSEstar2)
end
inc.totSqResidAll2 = norm(absresid(:,1).*absresid(:,1),1);
if maxResid > 0
    inc.numCloseAll2 = sum(absresid(:,1) <= maxResid);
    fprintf("  %d close points\n",inc.numCloseAll2)
end
inc.weights2 = w;
inc.w02 = w0;

% Update output solution if appropriate
if maxResid ~= 0
    if inc.numCloseAll2 >= inc.numCloseAll1
        outStep = 2;
        inc.weightsOut = inc.weights2;
        inc.w0Out = inc.w02;
        inc.residOut = inc.resid2;
        inc.totSqResidAllOut = inc.totSqResidAll2;
        inc.numCloseAllOut = inc.numCloseAll2;
        if mgood > 0
            inc.totSqResidTruOut = inc.totSqResidTru2;
            inc.numCloseTruOut = inc.numCloseTru2;
            inc.RMSETruOut = inc.RMSETru2;
            inc.MSETruOut = inc.MSETru2;
            inc.MAETruOut = inc.MAETru2;
            inc.TSEstarOut = inc.TSEstar2;
        end
    end
end

%--------------------------------------------------------------------------
% Reinstate any nonoutliers and regress again to get HP3

outliers = isoutlier(absresid);
% Make sure there are at least norig points left for the final fit
if m - sum(outliers) < norig
    outliers = zeros(m,1);
    sortedabsResid = sortrows([(1:m)',absresid],2,'descend');
    for i = 1:(m-norig)
        outliers(sortedabsResid(i,1),1) = 1;
    end
end

% Calculate outlier fractions at this stage (step 5)
if mgood > 0
    inc.HP3FracInOut = sum(outliers(1:mgood)/mgood);
    inc.HP3FracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end
   
B = zeros(m,norig);
y1 = zeros(m,1);
icount = 0;
for i=1:m
    if ~outliers(i,1)
        icount = icount + 1;
        B(icount,:) = Aorig(i,:);
        y1(icount,1) = y(i,1);
    end
end
B = B(1:icount,:);
y1 = y1(1:icount,1);
inc.m3 = icount;
fprintf("Reinstatement phase: %d points.\n",icount)

curr_time = toc;
fprintf("reinstated outliers %f\n", curr_time);

% Get the regression solution
A = [ones(icount,1),B];
beta = regress(y1,A);
w0 = beta(1,1);
w = beta(2:norig+1,1);
resid = w0 + Aorig*w - y;
absresid = abs(resid);

curr_time = toc;
fprintf("third regression %f\n", curr_time);

inc.resid3 = resid;
inc.totSqResidAll3 = norm(absresid(:,1).*absresid(:,1),1);
inc.weights3 = w;
inc.w03 = w0;

if mgood > 0
    inc.totSqResidTru3 = norm(absresid(1:mgood,1).*absresid(1:mgood,1),1);
    inc.RMSETru3 = sqrt(inc.totSqResidTru3);
    inc.MSETru3 = inc.totSqResidTru3/mgood;
    inc.MAETru3 = sum(abs(absresid(1:mgood,1)))/mgood;
    if maxResid > 0
        inc.numCloseTru3 = sum(absresid(1:mgood,1) <= maxResid);
    end
    % Calculate TSEstar
    sortedabsResid = sort(absresid);
    inc.TSEstar3 = norm(sortedabsResid(1:mgood,1).*sortedabsResid(1:mgood,1),1);
    fprintf("  TSEstar %f.\n",inc.TSEstar3)
end
if maxResid > 0
    inc.numCloseAll3 = sum(absresid(:,1) <= maxResid);
    fprintf("  %d close points.\n",inc.numCloseAll3)
end

if maxResid == 0
   % Just return the third regression
   outStep = 3;
   inc.weightsOut = inc.weights3;
   inc.w0Out = inc.w03;
   inc.residOut = inc.resid3;
   inc.totSqResidAllOut = inc.totSqResidAll3;
   if mgood > 0
       inc.totSqResidTruOut = inc.totSqResidTru3;
       inc.RMSETruOut = inc.RMSETru3;
       inc.MSETruOut = inc.MSETru3;
       inc.MAETruOut = inc.MAETru3;
       inc.TSEstarOut = inc.TSEstar3;
   end
end

%Update output solution if appropriate
if maxResid ~= 0
    if inc.numCloseAll3 >= inc.numCloseAllOut
        outStep = 3;
        inc.weightsOut = inc.weights3;
        inc.w0Out = inc.w03;
        inc.residOut = inc.resid3;
        inc.totSqResidAllOut = inc.totSqResidAll3;
        inc.numCloseAllOut = inc.numCloseAll3;
        if mgood > 0
            inc.totSqResidTruOut = inc.totSqResidTru3;
            inc.numCloseTruOut = inc.numCloseTru3;
            inc.RMSETruOut = inc.RMSETru3;
            inc.MSETruOut = inc.MSETru3;
            inc.MAETruOut = inc.MAETru3;
            inc.TSEstarOut = inc.TSEstar3;
        end
    end
end

%--------------------------------------------------------------------------
% Check phase: repeat the outlier reinstatement and removal and try again

outliers = isoutlier(absresid);
% Make sure there are at least norig points left for the final fit
if m - sum(outliers) < norig
    outliers = zeros(m,1);
    sortedabsResid = sortrows([(1:m)',absresid],2,'descend');
    for i = 1:(m-norig)
        outliers(sortedabsResid(i,1),1) = 1;
    end
end

% Calculate outlier fractions at this stage
if mgood > 0
    inc.HP4FracInOut = sum(outliers(1:mgood)/mgood);
    inc.HP4FracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end
   
B = zeros(m,norig);
y1 = zeros(m,1);
icount = 0;
for i=1:m
    if ~outliers(i,1)
        icount = icount + 1;
        B(icount,:) = Aorig(i,:);
        y1(icount,1) = y(i,1);
    end
end
curr_time = toc;
fprintf("reinstated outliers again %f\n", curr_time);

B = B(1:icount,:);
y1 = y1(1:icount,1);
inc.m4 = icount;
fprintf("Check phase: %d points.\n",icount)

% Get the regression solution
A = [ones(icount,1),B];
beta = regress(y1,A);
w0 = beta(1,1);
w = beta(2:norig+1,1);
resid = w0 + Aorig*w - y;
absresid = abs(resid);

curr_time = toc;
fprintf("regression again %f\n", curr_time);

inc.resid4 = resid;
inc.totSqResidAll4 = norm(absresid(:,1).*absresid(:,1),1);
inc.weights4 = w;
inc.w04 = w0;

if mgood > 0
    inc.totSqResidTru4 = norm(absresid(1:mgood,1).*absresid(1:mgood,1),1);
    inc.RMSETru4 = sqrt(inc.totSqResidTru4);
    inc.MSETru4 = inc.totSqResidTru4/mgood;
    inc.MAETru4 = sum(abs(absresid(1:mgood,1)))/mgood;
    if maxResid > 0
        inc.numCloseTru4 = sum(absresid(1:mgood,1) <= maxResid);
    end
    % Calculate TSEstar
    sortedabsResid = sort(absresid);
    inc.TSEstar4 = norm(sortedabsResid(1:mgood,1).*sortedabsResid(1:mgood,1),1);
    fprintf("  TSEstar %f.\n",inc.TSEstar4)
end
if maxResid > 0
    inc.numCloseAll4 = sum(absresid(:,1) <= maxResid);
    fprintf("  %d close points.\n",inc.numCloseAll4)
end

if maxResid == 0
   % Just return the 4th regression
   outStep = 4;
   inc.weightsOut = inc.weights4;
   inc.w0Out = inc.w04;
   inc.residOut = inc.resid4;
   inc.totSqResidAllOut = inc.totSqResidAll4;
   if mgood > 0
       inc.totSqResidTruOut = inc.totSqResidTru4;
       inc.RMSETruOut = inc.RMSETru4;
       inc.MSETruOut = inc.MSETru4;
       inc.MAETruOut = inc.MAETru4;
       inc.TSEstarOut = inc.TSEstar4;
   end
end

%Update output solution if appropriate
if maxResid ~= 0
    if inc.numCloseAll4 >= inc.numCloseAllOut
        outStep = 4;
        inc.weightsOut = inc.weights4;
        inc.w0Out = inc.w04;
        inc.residOut = inc.resid4;
        inc.totSqResidAllOut = inc.totSqResidAll4;
        inc.numCloseAllOut = inc.numCloseAll4;
        if mgood > 0
            inc.totSqResidTruOut = inc.totSqResidTru4;
            inc.numCloseTruOut = inc.numCloseTru4;
            inc.RMSETruOut = inc.RMSETru4;
            inc.MSETruOut = inc.MSETru4;
            inc.MAETruOut = inc.MAETru4;
            inc.TSEstarOut = inc.TSEstar4;
        end
    end
end

%--------------------------------------------------------------------------
% Collate the final output
curr_time = toc;
fprintf("calculating final output %f\n", curr_time);

inc.solTime = toc;
fprintf("Output solution from Step %d\n",outStep)
if mgood > 0
    fprintf("  TSEstar %f\n",inc.TSEstarOut);
    fprintf("  R %f\n",inc.TSEstarOut/inc.bnd2);
end
if maxResid ~= 0
    fprintf("  %d close points\n",inc.numCloseAllOut);
end

% Calculate the estimated number of inliers at output (inc.qout)
outliers = isoutlier(abs(inc.residOut));
inc.qout = m - sum(outliers);
% Calculate outlier fractions at final fit
if mgood > 0
    inc.FinalFracInOut = sum(outliers(1:mgood)/mgood);
    inc.FinalFracOutOut = sum(outliers(mgood+1:m)/(m-mgood));
end

return
end

%==========================================================================
% Analyzes the input data set to robustly identify outliers in regression 
% data sets.
%   Step 1 identifies outliers in the independent variables by looking
% for outliers in the original independent variable columns and then in
% the revised independent variable columns after a PCA transformation. A
% robust measure takes the largest value of the two and then looks for an
% abrupt change signalling the division between inliers and outliers.
%   Step 2 identifies outliers in the dependent variable. This can only be
% done by looking at the residuals relative to a regression fit. We use the
% residuals from the regression obtained over the set of inlier points 
% identified in Step 1.
%   Finally, we return the combined set of outliers identified in Steps 1
% and 2.
%   The robust measure in Step 1 scales the relative distance from the 
% median distance to the median for each point in each analysis above. Uses
% an offset in case the median (used as a denominator) is zero, or close.
%   If mgood > 0 then it returns measures on the accuracy of identifying
% outliers.
% INPUTS:
%   A: data matrix for independent variables
%   y: vector of dependent variable values
%   mgood: number of non-outliers, if known, in which case 
%          statistics are generated on how accurate the outlier
%          identification is. If mgood <= 0, or mgood >= m, then no such
%          statistics are generated. This is for testing purposes.
% OUTPUTS: these are all fields in outMeasure
%   .max: the maximum relative distance from the median, in medians, in any
%         of (i) through (iv) above, for each point
%   .TF: mx1 logical vector listing points identified as outliers
%   .count: total number of points identified as outliers
%   .q: the estimated number of inliers
%   If mgood > 0 and mgood < m:
%     .outTru: number of nonoutlier points identified as outliers
%     .outTruFrac: fraction of nonoutlier points identified as outliers
%     .outOut: number of outlier points identified as outliers
%     .outOutFrac: fraction of outlier points identified as outliers
function [outMeasure] = outFinderReg(A,y,mgood)
m = size(A,1);
n = size(A,2);
outMeasure.max = zeros(m,1) - Inf;
outMeasure.numZeroMedians = 0;

% columns of the independent variable matrix A
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

% principal component axes of the independent variable matrix A
[~,score,~] = pca(A);
for j=1:size(score,2)
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

% remove outliers in the independent variables
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
for i = istart:m
    if changes(i,1)
        cutoff = sortedMax(i,1);
        break
    end
end
indOutliers = outMeasure.max >= cutoff;
% eliminate those outliers and refit the PCA
B = zeros(m,n);
y1 = zeros(m,1);
row = zeros(m,1);
icount = 0;
for i=1:m
    if ~indOutliers(i,1)
        icount = icount + 1;
        B(icount,:) = A(i,:);
        y1(icount,1) = y(i,1);
        row(icount,1) = i;
    end
end
B = B(1:icount,:);
y1 = y1(1:icount,1);
row = row(1:icount,1);

% regress on the subset and see residuals over all points
A1 = [ones(icount,1),B];
beta = regress(y1,A1);
w0 = beta(1,1);
w = beta(2:n+1,1);
resid = w0 + A*w - y;
absresid = abs(resid);

% get the subset of residuals for the included points
absresidSub = absresid(indOutliers == 0,1);
% look for outliers or breakpoints
sortedARS = sort(absresidSub);
changes = ischange(sortedARS,'linear');
% Set initial defaults for the cutoff value. If there is no abrupt change
% then the default cutoff values are used.
if m/n <= 2
    % Can't start at m/2 or we don't leave enough points for PCA to run
    istart = n + 1;
    cutoff = sortedMax(istart,1);
else
    istart = ceil(m/2);
    if istart > size(sortedARS,1)
        cutoff = size(sortedARS,1);
    else
        cutoff = sortedARS(istart,1);
    end
end
for i = istart:size(absresidSub,1)
    if changes(i,1)
        cutoff = sortedARS(i,1);
        break
    end
end
depOutliers = absresidSub >= cutoff;

% Create one outlier list from the dependent and independent outliers
outliers = indOutliers;
for i=1:size(depOutliers,1)
    if depOutliers(i,1)
        outliers(row(i,1),1) = 1;
    end
end

outMeasure.TF = outliers;
outMeasure.count = sum(outMeasure.TF);
outMeasure.q = m - outMeasure.count;

% if mgood is nonzero, then calculate some accuracy measures
if (mgood > 0) && (mgood < m) 
    outMeasure.outTru = sum(outMeasure.TF(1:mgood,1));
    outMeasure.outTruFrac = outMeasure.outTru/mgood;
    outMeasure.outOut = sum(outMeasure.TF(mgood+1:m,1));
    outMeasure.outOutFrac = outMeasure.outOut/(m-mgood);
end

return
end
