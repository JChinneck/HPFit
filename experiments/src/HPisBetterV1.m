% Determines whether the hyperplane yielding edist2 is better than the
% hyperplane yielding edist1, or not.
% HP1 is the older hyperplane, HP2 is the newer one.
function [tf] = HPisBetterV1(edist1,edist2,maxDist,minFrac,feaTol)
m = size(edist1,1);
fprintf("  numEOutliers changes from %d to %d\n",sum(isoutlier(edist1)),sum(isoutlier(edist2)))

% Check the total number of points that are "close"
closeAll1 = sum(edist1(:,1) <= maxDist);
closeAll2 = sum(edist2(:,1) <= maxDist);
if closeAll2 > closeAll1
    fprintf("  CloseAll better: from %d to %d\n",closeAll1,closeAll2)
    tf = 1;
    return
end
if closeAll1 > closeAll2
    fprintf("  CloseAll worse: from %d to %d\n",closeAll1,closeAll2)
    tf = 0;
    return
end
if closeAll1 == closeAll2
    fprintf("  Closeall unchanged at %d.\n",closeAll1)
end

% There is a tie, so check the better ratios
numBetter = 0;
numWorse = 0;
for i=1:m
    if edist1(i,1) - edist2(i,1) > feaTol
        numBetter = numBetter + 1;
    end
    if edist2(i,1) - edist1(i,1) > feaTol
        numWorse = numWorse + 1;
    end
end
betterFrac = 0.0;
if numBetter > 0
    betterFrac = numBetter/(numBetter+numWorse);
end
fprintf("  numBetter %d numWorse %d better fraction %f.\n",numBetter,numWorse,betterFrac)

if (betterFrac >= minFrac) && (closeAll2 >= closeAll1)
    tf = 1;
    fprintf("  ***New hyperplane accepted.\n")
else
    tf = 0;
    fprintf("  ***New hyperplane rejected.\n")
end

return
end
