% Try reinstating points vs. some input edist vector and stillin vector.
% You can restate vs the current HP or vs. the incumbent HP as specified by
% vsWhich ('current' or 'incumbent'). You must set the value of vsWhich 
% before calling this routine. 

% See if some points are not strictly outliers distance-wise from the
% specified HP. If so, reinstate them and find the new hyperplane.

if strcmp(vsWhich,'incumbent')
    TF = isoutlier(errs.edist(:,3));
    stillin = errs.stillin(:,3);
else
    TF = isoutlier(edist);
end
numReinstated = 0;
B2 = zeros(m,norig);
stillin2 = stillin;
mnew = 0;
for i=1:m
    if stillin(i,1)
        mnew = mnew + 1;
        B2(mnew,:) = Aorig(i,:);
        continue
    end
    if ~TF(i,1)
        mnew = mnew + 1;
        B2(mnew,:) = Aorig(i,:);
        numReinstated = numReinstated + 1;
        stillin2(i,1) = 1;
    end
end
B2 = B2(1:mnew,:);
fprintf("  %d points reinstated. %d points total.\n",numReinstated,sum(stillin2))
if numReinstated > 0
    if dep_var == true
        beta = regress(B2(:,1),B2(:,2:n)); % regression with first variable as dependent
        beta_star = [-1.0; beta]; % coefficient on dependent variable is -1
        dist = abs(Aorig*beta_star); % get the absolute value of the difference between the dependent variable and the vertical projection on the regression hyperplane
        weights = beta; % for storing the coefficients; include the intercept; omit the response
        RHS = -beta(1); % just the intercept
        edist = dist;
    else
        % Find new hyperplane using PCA
        [weights,RHS] = getPCAHP(B2(:,2:n));
        % Calculate the point distances from the hyperplane
        dist = Aorig(:,2:n)*weights - RHS;
        % Calculate the Euclidean point distances from the hyperplane
        gradLen = norm(weights);
        edist = abs(dist/gradLen);
    end
    errs.totHPs = errs.totHPs + 1;
    
    % If the noBetter stopping condition is in use
    if strcmp(stopCondn,'noBetter')
        % See if the hyperplane is better. Previous best is always 3.
        if HPisBetterV1(errs.edist(:,3),edist(:,1),maxDist,minFrac,feaTol)
            gotBetter = 1;
            errs.bestHP = errs.totHPs;
            errs.weights(:,3) = weights;
            errs.RHS(1,3) = RHS;
            errs.stillin(:,3) = stillin2;
            errs.edist(:,3) = edist;
            errs.numEOutliers(1,3) = sum(isoutlier(edist));
            errs.totSqDistAll(1,3) = norm(errs.edist(:,3).*errs.edist(:,3),1);
            errs.numCloseAll(1,3) = sum(errs.edist(:,3) <= maxDist);
            errs.totRemoved(1,3) = m - sum(stillin2);
            sortedEdist = sort(edist,'ascend');
            errs.TSE(1,3) = sum(sortedEdist(1:TSEcount,1));
            if mgood > 0
                errs.totSqDistTru(1,3) = norm(errs.edist(1:mgood,3).*errs.edist(1:mgood,3),1);
                errs.numCloseTru(1,3) = sum(errs.edist(1:mgood,3) <= maxDist);
                errs.wrongRemoved(1,3) = mgood - sum(stillin2(1:mgood));
                errs.ratioTotSqDistTru(1,3) = 1.0;
                if errs.totSqDistTru(1,2) > 0.0
                    errs.ratioTotSqDistTru(1,3) = errs.totSqDistTru(1,3)/errs.totSqDistTru(1,2);
                end
                errs.diffCloseTru(1,3) = errs.numCloseTru(1,3) - errs.numCloseTru(1,2);
            end
            B = B2;
            stillin = stillin2;
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
        B = B2;
        stillin = stillin2;
    end
end
