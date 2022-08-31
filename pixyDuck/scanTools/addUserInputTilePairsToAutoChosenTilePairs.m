function tileIDpairsForAngle=addUserInputTilePairsToAutoChosenTilePairs(tileIDpairsForAngleUserInput,tileIDpairsForAngle,numTilePairsToTryToFindAngleFor)
% put user input pairs at top of list, remove duplicates if there are some, then trim total length to numTilePairsToTryToFindAngleFor or less
rng(0)

tileIDpairsForAngle=unique([tileIDpairsForAngleUserInput;tileIDpairsForAngle],'rows','stable');
numTilePairsToRemove=size(tileIDpairsForAngle,1)-numTilePairsToTryToFindAngleFor; % now with user input, may need to remove some tileIDpairs
numTilePairsToRemove(numTilePairsToRemove<0)=0;% if <0, make 0.

if numTilePairsToRemove>0
    iTilePairsToPotentiallyRemove=size(tileIDpairsForAngleUserInput,1)+1:size(tileIDpairsForAngle,1); % don't remove user input tile pairs
    
    if size(tileIDpairsForAngleUserInput,1)<numTilePairsToTryToFindAngleFor
        iTilePairsToKeep=[1:size(tileIDpairsForAngleUserInput,1),iTilePairsToPotentiallyRemove(sort(randperm(length(iTilePairsToPotentiallyRemove),length(iTilePairsToPotentiallyRemove)-numTilePairsToRemove)))];
    else
        iTilePairsToKeep=1:numTilePairsToTryToFindAngleFor;
    end
    tileIDpairsForAngle=tileIDpairsForAngle(iTilePairsToKeep,:);
end

end