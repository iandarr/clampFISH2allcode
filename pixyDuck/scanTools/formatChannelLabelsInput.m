function channelLabelsCellArray=formatChannelLabelsInput(Scan,channelLabelsInput)

numRounds=Scan.numRounds;

channelLabelsCellArray=repmat({{''}},numRounds,1);

if ~isempty(channelLabelsInput)
    if iscell(channelLabelsInput)
        if length(channelLabelsInput)~=numRounds
            error('user-input channelLabels should be a numRoundsx1 cell array of 1xnumChannels cell arrays')
        end
        
        for thisRound=1:Scan.numRounds
            channelLabelsThisRound=channelLabelsInput{thisRound};
            if ~iscell(channelLabelsThisRound)
                error('each of the elements of the user-input channelLabel cell array should itself be a cell array of size 1xnumChannels')
            elseif or(size(channelLabelsThisRound,2)~=Scan.Rounds(thisRound).numChannels,size(channelLabelsThisRound,1)~=1)
                chanNamesStr=join(Scan.Rounds(thisRound).channelNames,"', '");chanNamesStr=chanNamesStr{1};
                error('Element %i of user-input channelLabels should be 1x%i, one corresponding to each channel(%s) but it is %ix%i',thisRound,Scan.Rounds(thisRound).numChannels,chanNamesStr,size(channelLabelsThisRound,1),size(channelLabelsThisRound,2))
            else % store in Scan
                channelLabelsCellArray(thisRound)=channelLabelsThisRound;
            end
        end
    elseif istable(channelLabelsInput)
        Tc=channelLabelsInput;
        if ~all(ismember({'round','channel','channelLabel'},Tc.Properties.VariableNames))
            error("if channelLabels input is a table, it must have VariableNames 'round','channel',and 'channelLabel'")
        end
         for thisRound=1:numRounds
            allChannelNames=    Scan.Rounds(thisRound).channelNames;
            inputChannels=      Tc.channel(thisRound==Tc.round)';
            inputChannelLabels= Tc.channelLabel(thisRound==Tc.round)';
            %[a,idxReorder]=ismember(allChannelNames,inputChannels)
            channelLabelsThisRound=repmat({''},1,length(allChannelNames));
            for iChannel=1:length(allChannelNames)
                channel=allChannelNames{iChannel};
                idxInput=find(ismember(inputChannels,channel));
                if ~isempty(idxInput)
                    if ~isscalar(idxInput)
                        error('issue with input channelLabels table: should only have a single row for a given round and channel')
                    else
                        channelLabelsThisRound{iChannel}=inputChannelLabels{idxInput};
                    end
                end
            end
            channelLabelsCellArray{thisRound}=channelLabelsThisRound;
        end
        
    else
        error('unknown channelLabels input')
        
    end
else % input is isempty, return it empty
    for thisRound=1:numRounds
        numChannels=Scan.Rounds(thisRound).numChannels;
        channelLabelsCellArray{thisRound}=repmat({''},1,numChannels);
    end
end




end