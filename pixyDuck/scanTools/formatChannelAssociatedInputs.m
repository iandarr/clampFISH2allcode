function channelPropTable=formatChannelAssociatedInputs(Scan,channelLabelsInput,channelLabelsExactMatchOnly)
warning('off','MATLAB:table:RowsAddedExistingVars')
numRounds=Scan.numRounds;
channelPropTable=table();%array2table('Size',[numRounds,1],'VariableNames','channelLabel');

channelLabelsCellArray=repmat({{''}},numRounds,1);

if ~isempty(channelLabelsInput)
    if iscell(channelLabelsInput)
        % should just be single columnn cell array with channelLabel in it
        if length(channelLabelsInput)~=numRounds
            error('user-input channelLabels should be a numRoundsx1 cell array of 1xnumChannels cell arrays')
        end
        
        for thisRound=1:Scan.numRounds
            channelPropValsThisRound=channelLabelsInput{thisRound};
            if ~iscell(channelPropValsThisRound)
                error('each of the elements of the user-input channelLabels cell array should itself be a cell array of size 1xnumChannels')
            elseif or(size(channelPropValsThisRound,2)~=Scan.Rounds(thisRound).numChannels,size(channelPropValsThisRound,1)~=1)
                chanNamesStr=join(Scan.Rounds(thisRound).channelNames,"', '");chanNamesStr=chanNamesStr{1};
                error('Element %i of user-input channelLabels should be 1x%i, one corresponding to each channel(%s) but it is %ix%i',thisRound,Scan.Rounds(thisRound).numChannels,chanNamesStr,size(channelPropValsThisRound,1),size(channelPropValsThisRound,2))
            else % store in Scan
                channelLabelsCellArray(thisRound)=channelPropValsThisRound;
            end
        end
        channelPropTable.channelLabels=channelLabelsCellArray;
        
    elseif istable(channelLabelsInput)
        Tc=channelLabelsInput;
        if ~all(ismember({'round','channel','channelLabels'},Tc.Properties.VariableNames))
            error("if channelLabels input is a table, it must have VariableNames 'round','channel',and 'channelLabels'")
        end
        
        channelProps=Tc.Properties.VariableNames(~ismember(Tc.Properties.VariableNames,{'round','channel'}));
        
        for thisRound=1:numRounds
            allChannelNames=    Scan.Rounds(thisRound).channelNames;
            inputChannels=      Tc.channel(thisRound==Tc.round)';
            
            for thisProp=1:length(channelProps)
                prop=channelProps{thisProp};
                thisPropInput= Tc.(prop)(thisRound==Tc.round)';
                %[a,idxReorder]=ismember(allChannelNames,inputChannels)
                %propClass=class(Tc.(inputProp));
                
                
                if isequal(class(Tc.(prop)),'cell')
                    channelPropValsThisRound=repmat({''},1,length(allChannelNames));
                    propIsNumeric=false;
                elseif isnumeric(Tc.(prop))
                    channelPropValsThisRound=nan(1,length(allChannelNames));
                    propIsNumeric=true;
                else
                    error('channelLabelsInput columns must be numeric or cell')
                end
                
                for iChannel=1:length(allChannelNames)
                    channel=allChannelNames{iChannel};
                    if channelLabelsExactMatchOnly
                        idxInput=ismember(inputChannels,channel); % exact match
                        if sum(idxInput)>1
                            error('issue with input channelLabels table: should only have a single row for a given round and channel')
                        end
                    else
                        idxInput=arrayfun(@(x) startsWith(channel,x),inputChannels);%startsWith
                        if sum(idxInput)>1
                            error('issue with input channelLabels table using channelLabelsExactMatchOnly==false: should only have a single row applicable for a given round and channel')
                        end
                    end
                    
                    if sum(idxInput)==1
                        if propIsNumeric
                            channelPropValsThisRound(iChannel)=thisPropInput(idxInput);
                        else
                            channelPropValsThisRound{iChannel}=thisPropInput{idxInput};
                        end
                    end
                end
                
                
                channelPropTable.(prop)(thisRound)={channelPropValsThisRound};
                %channelLabelsCellArray{thisRound}=channelLabelsThisRound;
                
            end
            
            
        end
        
    else
        error('unknown channelLabels input')
        
    end
else % input is isempty, return it empty
    for thisRound=1:numRounds
        numChannels=Scan.Rounds(thisRound).numChannels;
        channelLabelsCellArray{thisRound}=repmat({''},1,numChannels);
    end
    channelPropTable.channelLabels=channelLabelsCellArray;
end



warning('on','MATLAB:table:RowsAddedExistingVars')
end