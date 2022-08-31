function checkControlPointInitialGuessesTable(T,numRounds,refRound)
% check input table controlPointInitialGuessesTable
%
% called T here for clarity

if isempty(T)
    return
    % empty table is okay
end

if or(~istable(T),~all(ismember({'fixedRound','moveRound','fixedX','fixedY','moveX','moveY'},T.Properties.VariableNames)))
    error('input controlPointInitialGuessesTable must be a table with VariableNames fixedRound moveRound fixedX fixedY moveX moveY')
end

if ~all([isnumeric(T.fixedRound),isnumeric(T.moveRound),isnumeric(T.fixedX),isnumeric(T.fixedY),isnumeric(T.moveX),isnumeric(T.moveY)])
    error('input controlPointInitialGuessesTable has at least one non-numeric type out of variables: fixedRound moveRound fixedX fixedY moveX moveY')
end

if ~all(ismember([T.fixedRound;T.moveRound],1:numRounds))
    error('input controlPointInitialGuessesTable references rounds that are not in the range of [1:%i]',numRounds)
end

if any(T.moveRound==refRound)
    error('input controlPointInitialGuessesTable variable moveRound cannot include the referenceRound (%i), since it is by definition not being moved',refRound)
end

if any(T.fixedRound==T.moveRound)
    error('input controlPointInitialGuessesTable issue: a single row must have a different fixedRound and moveRound')
end

moveRoundsIncluded=unique(T.moveRound);
for i=1:length(moveRoundsIncluded)
    thisMoveRound=moveRoundsIncluded(i);
    idxsThisMoveRound=find(T.moveRound==thisMoveRound);
    if ~all(T.fixedRound(idxsThisMoveRound)==T.fixedRound(idxsThisMoveRound(1)))
        error('input controlPointInitialGuessesTable issue: moveRound %i has more than one fixedRound associated with it. It can only have one',thisMoveRound)
    end
end

end