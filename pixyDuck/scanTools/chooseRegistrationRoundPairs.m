function pairsOut=chooseRegistrationRoundPairs(numRounds,refRound,pairsInput,varargin)
% pairsOut=chooseRegistrationRoundPairs(numRounds,refRound,pairsInput)
% pairsOut=chooseRegistrationRoundPairs(numRounds,refRound,pairsInput,controlPointInitialGuessesTable)
%
% numRounds=Scan.numRounds;
% refRound=Scan.referenceRound;

if nargin>3
    if nargin==4
        controlPointInitialGuessesTable=varargin{1};
        checkControlPointInitialGuessesTable(controlPointInitialGuessesTable,numRounds,refRound) %check over this input
    else
        error('cannot interpret last input argument')
    end
    adjustPairsForControlPointInitialGuesses=true;
else
    adjustPairsForControlPointInitialGuesses=false;
end

if isempty(pairsInput) % make default registrationRoundPairs
    
    pairs=makeDefaultRefistrationRoundPairs(refRound,numRounds);
    
else % make the pairsOut including the user registrationRoundPairs (pairsInput)
    if ~isnumeric(pairsInput)
        error('input registrationRoundPairs must be a numeric matrix with 2 columns, with first column specifying the fixed rounds and the second column specifying the round to move')
    elseif ~(size(pairsInput,2)==2)
        error('input registrationRoundPairs must have 2 columns')
    end
    
    numPairsInput=size(pairsInput,1);
    % check over pairsInput
    if any(pairsInput(:,2)==refRound)
        error('second column of registrationRoundPairs cannot include the reference round (%i)',refRound)
    end
    
    if length(pairsInput(:,2))~=length(unique(pairsInput(:,2)))
        error('second column of registrationRoundPairs should only include each round once')
    end
    if any(pairsInput(:,1)==pairsInput(:,2))
        error('registrationRoundPairs input cannot have the same round in columns 1 and 2')
    end
    
    if ~all(ismember([pairsInput(:);refRound],1:numRounds))
        error('registrationRoundPairs input and refRound must all be integers 1 or greater')
    end
    
    if rem(numRounds,1)~=0
        error('numRounds input must be integer')
    end
    
    % then make registrationRoundPairs with this input
    pairs=makeDefaultRefistrationRoundPairs(refRound,numRounds);
    pairs=pairs(~(ismember(pairs(:,2),pairsInput(:,2))),:); % take out default rows related to rounds in inputs
    pairsInputRemaining=pairsInput; %initialize
    for i=1:size(pairsInput,1)
        roundsWithRegistration=unique([refRound;pairs(:,2)]);
        idxPairToInclude=find(ismember(pairsInputRemaining(:,1),roundsWithRegistration),1,'first');
        if isempty(idxPairToInclude)
            error('could not find a complete registrationRoundPairs matrix that satisfied the input provided. Check that two of the input rows provided arent referencing eachother in a circular manner')
        end
        pairs=[pairs;pairsInputRemaining(idxPairToInclude,:)];
        pairsInputRemaining=pairsInputRemaining([1:size(pairsInputRemaining,1)]~=idxPairToInclude,:); %update
    end
    %pairsOut=pairs;
end

%% If there is a controlPointInitialGuessesTable, will try to re-arrange rows of pairsOut to accomodate these

%pairs
if adjustPairsForControlPointInitialGuesses
    
    
    numCpPairs=height(controlPointInitialGuessesTable);
    %controlPointInitialGuessesTable
    %pairs
    for i=1:numCpPairs
        % fixedR needs to be found lower down in col2 of pairs variable
        fixedR=controlPointInitialGuessesTable.fixedRound(i);
        moveR=controlPointInitialGuessesTable.moveRound(i);
        if fixedR==refRound
            % no problem with this, leave pairs as is
        else
            rowFixedR=find(pairs(:,2)==fixedR);
            rowMoveR=find(pairs(:,2)==moveR);
            if rowFixedR<rowMoveR
                %3          2
                % no problem with this,
            else
                %fprintf('  need to reorder pairs matrix\n')
                % try to take row rowFixedR and put it at row rowMoveR, shift everything else down
                newRowOrder=[1:rowMoveR-1,rowFixedR,rowMoveR:size(pairs,1)];
                %newRowOrder=newRowOrder(~find(newRowOrder==rowFixedR,1,'last'));
                idxToRemove=find(newRowOrder==rowFixedR,1,'last');
                newRowOrder=newRowOrder(~ismember(1:length(newRowOrder),idxToRemove)); %remove duplicate
                
                pairsCandidate=pairs(newRowOrder,:);
                
                % check if pairs matrix is okay
                nonRefRounds=[1:refRound-1,refRound+1:numRounds];
                for ii=1:length(nonRefRounds)
                    nonRefRound=nonRefRounds(ii);
                    firstRowAsFixedR=find(nonRefRound==pairsCandidate(:,1));
                    firstRowAsMoveR=find(nonRefRound==pairsCandidate(:,2));
                    if any(firstRowAsFixedR<=firstRowAsMoveR)
                        error('could not find a compatible registrationRoundPairs ordering with controlPointInitialGuessesTable row %i (fixedRound = %i, moveRound= %i). Check that these two inputs do not reference rounds in a circular manner',i,fixedR,moveR)
                    end
                    
                end
                % if we survived that it's okay
                pairs=pairsCandidate;
                
            end
        end
        
    end
    
    pairsOut=pairs;
else
    pairsOut=pairs;
end
end

function defaultPairs=makeDefaultRefistrationRoundPairs(refRound,numRounds)
temp=[1:numRounds]'; temp=temp(~(temp==refRound));
defaultPairs=[repmat(refRound,numRounds-1,1),temp];

end