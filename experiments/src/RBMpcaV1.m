% RBMpca uses PCA for hyperplane placement and standard outlier detection
% methods (3 median absolute distances criterion).
%    Purpose: to find the best fitting hyperplane to a set of input data,
% where best fitting means more points are closer to the hyperplane than
% competing alternatives. Note that this cardinality criterion differs
% quite a bit from measures such as total squared distance to the
% hyperplane, which can be fooled by outliers.
%    Specifically designed to be robust to outliers, especially clusters of
% outliers. The main idea to iteratively remove outliers, then place and
% evaluate intermediate hyperplanes until stopping conditions are met.
%    There is no assumption that one column is the target. It is looking
% for the best fit over all points and all variables.

% Main versions and options:
% Stopping condition:
%   - noCand: only when the outlier removal method detects no more outliers
%   - noBetter: when the "better" criterion decides that current HP is not
%       better than last HP.
% maxDist:
%   - part of the "better" criterion (see Notes below). This can be
%     prespecified, or found automatically. If negative, it specifies the
%     percentile distance from HP1 to define maxDist. If positive it
%     specifes the maxDist directly.
% Reinstate (0/1):
%   - whether or not to reinstate points that are not identified as
%     outliers relative to latest HP using Matlab isoutlier().
% mgood:
%   - for testing purposes. If there is a known hyperplane, then the first
%     mgood points are associated with it. If mgood <= 0 then there
%     is no a priori known hyperplane. If the mgood points are known, then
%     statistics about those points can be generated. Useful for testing.
% dep_var:
%   - true: a dependent variable is specified, as in ordinary 
%           regression.  A residual is measured as the absolute
%           difference between the response value and the
%           vertical projection of the point on the fitted hyperplane
%   - false: no dependent variable is specified. PCA is used as a 
%            subroutine to fit hyperplanes.  If insufficient 
%            data are available, an elastic LP is solved.
%            A residual is measured using the distance to the 
%            orthogonal projection.

% NOTES:
% - The "better" criterion works as follows. HP B is better than HP A if it
%   has more points that are relatively closer to HP B than HP A, and at
%   least as many points that are within maxDist Euclidean distance as for
%   HP A.
% - "traditional" method is a single outlier removal step followed by a
%   single HP placement step.
% - "iterative traditional" method is an iterative removal of outliers
%   followed by HP placement until no more outliers are detected.
% - the Relative Better Method (RBM) is iterative removal of outliers, with
%   reinstatement, nobetter exit, and automatic selection of maxDist
% - the data file is size mtot rows by norig columns.

% INPUTS:
%   parameters mentioned above
%   Aorig: raw data table
%   mgood: number of good rows of data. After that they are outliers
%   TSEcount: if a trimmed sum of absolute errors report is desired, set
%      this paramter to the required count. If <= 0 it is ignored.
%   RBMparam: the parameters controlling the solves
%     - feaTol: feasibility tolerance or equality tolerance
%         (default 1.0e-6)
%     - minFrac: new HP must have at least this fraction of changed
%         distances closer than farther away (default 0.5)
%     - maxRemoveFrac: can remove no more than this fraction of points
%        before abandoning any more outlier removals (default 0.75)
%
% OUTPUTS:
%   weights: the hyperplane weights.
%   RHS: the RHS constant.
%   numRemoved: total number of outliers removed in final hyperplane
%   solTime: solution time (s)
%   nHPs: number of HPs placed
%   errs: error measures for retained pts (ret) and all pts (all)
%     retTotAbsDist, allTotAbsDist: sum of absolute distances
%     retTotSqDist, allTotSqDist: sum of squared distances
%     retMedAbsDist, allMedAbsDist: median absolute distance
%     retMaxAbsDist, allMaxAbsDist: maximum absolute distance
%     outInitDist, outFinalDist: average distance from initial and final
%       hyperplanes for points deemed outliers in final hyperplane
%     betterFrac: fraction of points closer to final hyperplane than
%       initial hyperplane
% if mgood > 0 then these values are also reported:
%   - numWrong: total of non-outliers removed pts + true outliers not removed
%   - For the true points: errs.
%       - truTotAbsDist: sum of absolute distances
%       - truTotSqDist: sum of squared distances
%       - truMedAbsDist: median absolute distance
%       - truMaxAbsDist: maximum absolute distance

function [errs,solTime,nHPs] = RBMpcaV1(Aorig,mgood,TSEcount,RBMparam,dep_var)


% Get control parameter values
feaTol = RBMparam.feaTol;
if feaTol <= 0
    feaTol = 1.0e-6;
end
maxDist = RBMparam.maxDist;
minFrac = RBMparam.minFrac;
if minFrac <= 0
    minFrac = 0.5;
end
maxRemoveFrac = RBMparam.maxRemoveFrac;
if maxRemoveFrac <= 0
    maxRemoveFrac = 0.5;
end
stopCondn = RBMparam.stopCondn;
reinstate = RBMparam.reinstate;

tic;

% Get data table dimensions
m = size(Aorig,1);
errs.mtot = m;
errs.lts = 0.0;
norig = size(Aorig,2);
n = norig;
fprintf("mgood %d mtot %d n %d\n",mgood,m,norig)

% Initialize some data.
% 1: initial HP, 2: HP after single removal step, 3: final HP
errs.weights = zeros(norig-1,3); % minus 1 because of column of 1s
errs.RHS = zeros(1,3);
errs.stillin = [ones(m,1),zeros(m,2)];
errs.edist = zeros(m,3);
errs.numEOutliers = zeros(1,3);
errs.totSqDistTru = zeros(1,3) - 1;
errs.numCloseTru = zeros(1,3) - 1;
errs.totSqDistAll = zeros(1,3);
errs.numCloseAll = zeros(1,3);
errs.totRemoved = zeros(1,3);
errs.wrongRemoved = zeros(1,3);
errs.fracCloseTru = zeros(1,3);
errs.fracCloseAll = zeros(1,3);
errs.success = zeros(1,3); % -1:failure, 0:neutral, 1:success
errs.TSE = zeros(1,3) + Inf; % trimmed sum of edists
errs.gamma = zeros(1,3) + Inf; %qth residual

% PART 1: Initial Hyperplane, no removals ---------------------------------

% get the PCA solution for the original dataset
if dep_var == true
    fprintf("Using regression\n")
    beta = regress(Aorig(:,1), Aorig(:,2:n)); %regression with first variable as dependent
    beta_star = [-1.0; beta]; % coefficient on dependent variable is -1
    dist = abs(Aorig*beta_star); % get the absolute value of the difference between the dependent variable and the vertical projection on the regression hyperplane
    errs.weights(:,1) = beta; % for storing the coefficients, including the intercept, omit the response
    errs.RHS(:,1) = -beta(1); % just the intercept
    edist=dist;
else
    [errs.weights(:,1),errs.RHS(1,1)] = getPCAHP(Aorig(:,2:n));
    % Calculate the point distances from the hyperplane
    dist = Aorig(:,2:n)*errs.weights(:,1) - errs.RHS(1,1);
    % Calculate the Euclidean point distances from the hyperplane
    gradLen = norm(errs.weights(:,1));
    edist = abs(dist/gradLen);
end
errs.numEOutliers(1,1) = sum(isoutlier(edist));
% Find maxDist automatically if not prespecified
if maxDist <= 0
    maxDist = prctile(edist,-maxDist);
    fprintf("maxDist automatically selected as %f\n",maxDist)
end
errs.maxDist = maxDist; % record the final value of maxDist
errs.edist(:,1) = edist;
if mgood > 0
    errs.totSqDistTru(1,1) = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
    errs.numCloseTru(1,1) = sum(edist(1:mgood,1) <= maxDist);
end
errs.totSqDistAll(1,1) = norm(edist(:,1).*edist(:,1),1);
errs.numCloseAll(1,1) = sum(edist(:,1) <= maxDist);
fprintf("HP1 tot sq dist tru %f. Close: all %d tru %d. numEOutliers %d\n", ...
    errs.totSqDistTru(1,1),errs.numCloseAll(1,1),errs.numCloseTru(1,1), ...
    errs.numEOutliers(1,1))
sortedEdist = sort(edist,'ascend');
fprintf("TSEcount %d\n", TSEcount)
disp(TSEcount)
if TSEcount > 0
    errs.TSE(1,1) = sum(sortedEdist(1:TSEcount,1));
    errs.gamma(1,1) = sortedEdist(TSEcount);
    errs.lts = sum(sortedEdist(1:mgood,1).*sortedEdist(1:mgood,1));
end

% PART 2: Iterate to place improving hyperplanes --------------------------

% initialize HP2 (the traditional result) = HP1
errs.edist(:,2) = errs.edist(:,1);
errs.numEOutliers(1,2) = errs.numEOutliers(1,1);
errs.weights(:,2) = errs.weights(:,1);
errs.RHS(1,2) = errs.RHS(1,1);
errs.stillin(:,2) = errs.stillin(:,1);
errs.totSqDistTru(1,2) = errs.totSqDistTru(1,1) ;
errs.numCloseTru(1,2) = errs.numCloseTru(1,1);
errs.totSqDistAll(1,2) = errs.totSqDistAll(1,1);
errs.numCloseAll(1,2) = errs.numCloseAll(1,1);
errs.totRemoved(1,2) = errs.totRemoved(1,1);
errs.wrongRemoved(1,2) = errs.wrongRemoved(1,1);
errs.ratioTotSqDistTru(1,2) = 1.0;
errs.diffCloseTru(1,2) = 0;
errs.TSE(1,2) = errs.TSE(1,1);
errs.gamma(1,2) = errs.gamma(1,1);

% Initialize HP3 (output hyperplane) = HP1
errs.bestHP = 1;
errs.edist(:,3) = errs.edist(:,1);
errs.numEOutliers(1,3) = errs.numEOutliers(1,1);
errs.weights(:,3) = errs.weights(:,1);
errs.RHS(1,3) = errs.RHS(1,1);
errs.stillin(:,3) = errs.stillin(:,1);
errs.totSqDistTru(1,3) = errs.totSqDistTru(1,1) ;
errs.numCloseTru(1,3) = errs.numCloseTru(1,1);
errs.totSqDistAll(1,3) = errs.totSqDistAll(1,1);
errs.numCloseAll(1,3) = errs.numCloseAll(1,1);
errs.totRemoved(1,3) = errs.totRemoved(1,1);
errs.wrongRemoved(1,3) = errs.wrongRemoved(1,1);
errs.ratioTotSqDistTru(1,3) = 1.0;
errs.diffCloseTru(1,3) = 0;
errs.TSE(1,3) = errs.TSE(1,1);
errs.gamma(1,3) = errs.gamma(1,1);
errs.totHPs = 1;

% Start the loop of removals and improving hyperplanes
stillin = ones(m,1);
laststillin = zeros(m,1);
B = Aorig;
bestABMcount = 0;
ABMticker = 0;
for itn = 1:m
    
    if stillin == laststillin
        fprintf("Dataset unchanged since last iteration: exiting.\n")
        break
    end
    laststillin = stillin;
    
    if errs.numCloseAll(1,3) > bestABMcount
        bestABMcount = errs.numCloseAll(1,3);
        ABMticker = 0;
        fprintf("  ABM increases to %d at itn %d.\n",bestABMcount,itn)
    else
        ABMticker = ABMticker + 1;
        if ABMticker > 10
            % ABM count hasn't improved in 10 iterations, so break
            fprintf("  No change in ABM for 10 iterations: exiting.\n")
            break
        end
    end
    
    gotBetter = 0;
    if strcmp(stopCondn,'noCand')
        % If stopping only if there are no candidates, then always assume
        % that you found a better solution in each iteration until there
        % are no candidates
        gotBetter = 1;
    end
    
    % Track which rows are in B
    trackRows = find(stillin);
    
    % Remove outliers
    [B,removed] = rmoutliers(B);
    numRemoved = sum(removed);
    if numRemoved > 0
        for i=1:size(removed,1)
            if removed(i,1)
                stillin(trackRows(i,1),1) = 0;
            end
        end
    end
    fprintf("  %d outliers removed by stdOutliers. %d total points.\n",numRemoved,sum(stillin))
    
    % If outliers were removed, then find a new HP.
    % The noCand stopping condition assumes it always finds a better new
    % HP until there are no more candidates. This can be overridden by the
    % noBetter stopping condition.
    
    if numRemoved > 0
        if dep_var == true
            fprintf("Using regression\n")
            beta = regress(B(:,1), B(:,2:n)); % regression with first variable as dependent
            beta_star = [-1.0; beta]; % coefficient on dependent variable is -1; including the intercept
            dist = abs(Aorig*beta_star); % get the absolute value of the difference between the dependent variable and the vertical projection on the regression hyperplane
            weights = beta;% for storing the coefficients ; including the intercept; omit the response
            RHS = -beta(1); % just the intercept
            edist = dist;
        else
            % Find new hyperplane via PCA
            [weights,RHS] = getPCAHP(B(:,2:n));
            % Calculate the point distances from the hyperplane
            dist = Aorig(:,2:n)*weights - RHS;
            % Calculate the Euclidean point distances from the hyperplane
            gradLen = norm(weights);
            edist = abs(dist/gradLen);
        end
        %beta_star = beta_star(:,1);
        errs.totHPs = errs.totHPs + 1;
        
        if itn == 1
            % Must update traditional result (HP2)
            errs.weights(:,2) = weights;
            errs.RHS(1,2) = RHS;
            errs.stillin(:,2) = stillin;
            size(edist,1)
            errs.edist(:,2) = edist;
            errs.numEOutliers(1,2) = sum(isoutlier(edist));
            errs.totSqDistAll(1,2) = norm(edist(:,1).*edist(:,1),1);
            errs.numCloseAll(1,2) = sum(edist(:,1) <= maxDist);
            errs.totRemoved(1,2) = m - sum(stillin);
            sortedEdist = sort(edist,'ascend');
            errs.TSE(1,2) = sum(sortedEdist(1:TSEcount,1));
            errs.gamma(1,2) = sortedEdist(TSEcount);
            errs.lts = sum(sortedEdist(1:mgood,1).*sortedEdist(1:mgood,1));
            if mgood > 0
                errs.totSqDistTru(1,2) = norm(edist(1:mgood,1).*edist(1:mgood,1),1);
                errs.numCloseTru(1,2) = sum(edist(1:mgood,1) <= maxDist);
                errs.wrongRemoved(1,2) = mgood - sum(stillin(1:mgood));
                errs.ratioTotSqDistTru(1,2) = 1.0;
                if errs.totSqDistTru(1,1) > 0.0
                    errs.ratioTotSqDistTru(1,2) = errs.totSqDistTru(1,2)/errs.totSqDistTru(1,1);
                end
                errs.diffCloseTru(1,2) = errs.numCloseTru(1,2) - errs.numCloseTru(1,1);
            end
            fprintf("HPtrad tot sq dist tru %f. Close: all %d tru %d. numEOutliers %d\n", ...
                errs.totSqDistTru(1,2),errs.numCloseAll(1,2),errs.numCloseTru(1,2), ...
                errs.numEOutliers(1,2))
        end
        
        % if using the noBetter stopping condition, check the new HP
        if strcmp(stopCondn,'noBetter')
            % See if the hyperplane is better. Previous best is always 3.
            if HPisBetterV1(errs.edist(:,3),edist(:,1),maxDist,minFrac,feaTol)
                gotBetter = 1;
                errs.bestHP = errs.totHPs;
                errs.weights(:,3) = weights;
                errs.RHS(1,3) = RHS;
                errs.stillin(:,3) = stillin;
                errs.edist(:,3) = edist;
                errs.numEoutliers(1,3) = sum(isoutlier(edist));
                errs.totSqDistAll(1,3) = norm(errs.edist(:,3).*errs.edist(:,3),1);
                errs.numCloseAll(1,3) = sum(errs.edist(:,3) <= maxDist);
                errs.totRemoved(1,3) = m - sum(stillin);
                sortedEdist = sort(edist,'ascend');
                errs.TSE(1,3) = sum(sortedEdist(1:TSEcount,1));
                errs.lts = sum(sortedEdist(1:mgood,1).*sortedEdist(1:mgood,1));
                errs.gamma(1,3) = sortedEdist(TSEcount);
                if mgood > 0
                    errs.totSqDistTru(1,3) = norm(errs.edist(1:mgood,3).*errs.edist(1:mgood,3),1);
                    errs.numCloseTru(1,3) = sum(errs.edist(1:mgood,3) <= maxDist);
                    errs.wrongRemoved(1,3) = mgood - sum(stillin(1:mgood));
                    errs.ratioTotSqDistTru(1,3) = 1.0;
                    if errs.totSqDistTru(1,2) > 0.0
                        errs.ratioTotSqDistTru(1,3) = errs.totSqDistTru(1,3)/errs.totSqDistTru(1,2);
                    end
                    errs.diffCloseTru(1,3) = errs.numCloseTru(1,3) - errs.numCloseTru(1,2);
                end
            end
        else
            errs.bestHP = errs.totHPs;
            errs.weights(:,3) = weights;
            errs.RHS(1,3) = RHS;
            errs.stillin(:,3) = stillin;
            errs.edist(:,3) = edist;
            errs.numEoutliers(1,3) = sum(isoutlier(edist));
            errs.totSqDistAll(1,3) = norm(errs.edist(:,3).*errs.edist(:,3),1);
            errs.numCloseAll(1,3) = sum(errs.edist(:,3) <= maxDist);
            errs.totRemoved(1,3) = m - sum(stillin);
            if mgood > 0
                errs.totSqDistTru(1,3) = norm(errs.edist(1:mgood,3).*errs.edist(1:mgood,3),1);
                errs.numCloseTru(1,3) = sum(errs.edist(1:mgood,3) <= maxDist);
                errs.wrongRemoved(1,3) = mgood - sum(stillin(1:mgood));
                errs.ratioTotSqDistTru(1,3) = 1.0;
                if errs.totSqDistTru(1,2) > 0.0
                    errs.ratioTotSqDistTru(1,3) = errs.totSqDistTru(1,3)/errs.totSqDistTru(1,2);
                end
                errs.diffCloseTru(1,3) = errs.numCloseTru(1,3) - errs.numCloseTru(1,2);
            end
        end
    else
        % No point proceeding since no points removed
        fprintf("  No points removed. Exiting\n")
        break
    end   % on numRemoved
    
    % If allowing reinstatement check:
    if reinstate > 0
        vsWhich = 'current';
        tryReinstatingV1;
    end
    
    if gotBetter > 0
        fprintf("HP%d tot sq dist tru %f. Close: all %d tru %d. numEOutliers %d\n", ...
            errs.totHPs, errs.totSqDistTru(1,3),errs.numCloseAll(1,3),errs.numCloseTru(1,3), ...
            errs.numEOutliers(1,3))
    else
        fprintf("  No better solution found at iteration %d. Reinstating vs. incumbent.\n",itn)
        vsWhich = 'incumbent';
        tryReinstatingV1;
        if gotBetter > 0
            fprintf("HP%d tot sq dist tru %f. Close: all %d tru %d. numEOutliers %d\n", ...
                errs.totHPs, errs.totSqDistTru(1,3),errs.numCloseAll(1,3),errs.numCloseTru(1,3), ...
                errs.numEOutliers(1,3))
        else
            fprintf("  No better solution from incumbent restatement. Exiting.\n")
            break
        end
    end
    
    if sum(stillin)/m < 1.0 - maxRemoveFrac
        fprintf("  attempt to remove too many points. Bail out\n")
        break;
    end
    
end % of itn loop

solTime = toc;
nHPs = errs.totHPs;

fprintf("Output hyperplane: %d of %d\n",errs.bestHP,errs.totHPs)
errs.fracCloseAll(1,3) = errs.numCloseAll(1,3)/m;
errs.fracCloseAll(1,2) = errs.numCloseAll(1,2)/m;
errs.fracCloseAll(1,1) = errs.numCloseAll(1,1)/m;
fprintf("Total squared errors for all points. Initial %f. Final %f.\n", ...
    errs.totSqDistAll(1,1),errs.totSqDistAll(1,3))
if mgood > 0
    fprintf("Total squared errors for true points. Initial %f. Final %f.\n", ...
        errs.totSqDistTru(1,1),errs.totSqDistTru(1,3))
    errs.fracCloseTru(1,3) = errs.numCloseTru(1,3)/mgood;
    errs.fracCloseTru(1,2) = errs.numCloseTru(1,2)/mgood;
    errs.fracCloseTru(1,1) = errs.numCloseTru(1,1)/mgood;
    for ii=1:3
        fprintf("(%d) Fraction close all: %f tru: %f\n",ii, ...
            errs.fracCloseAll(1,ii),errs.fracCloseTru(1,ii))
    end
    errs.success = 0;
    if errs.totSqDistTru(1,3) < min(errs.totSqDistTru(1,1),errs.totSqDistTru(1,2)) - feaTol
        fprintf("SUCCESS!\n")
        errs.success = 1;
        return
    end
    if errs.totSqDistTru(1,3) > min(errs.totSqDistTru(1,1),errs.totSqDistTru(1,2)) + feaTol
        fprintf("FAILURE.\n")
        errs.success = -1;
    end
else
    fprintf("Fraction close all: %f\n",errs.fracCloseAll)
end

return

end
