function sTiffNames=tifNameFormatInput2TifNameStruct(Scan,tifNameFormat)
numRounds=Scan.numRounds;
sTiffNames=struct();
for thisRound=1:numRounds
    channelNames=Scan.Rounds(thisRound).channelNames;
    numChannels=length(channelNames);
    namedScanFields=fields(Scan.Rounds(thisRound));
    
    tifNamesThisRound=repmat({''},1,numChannels);
    %tifNameFormat={'R','round','_','channelName','_','channelLabel'};
    for iFormatSpec=1:length(tifNameFormat)
        formatSpec=tifNameFormat{iFormatSpec};
        
        if ~ischar(formatSpec)
            error('makeStitches: each cell of tifNameFormat should be a character vector')
        end
        
        if strcmp(formatSpec,'round')
            tifNamesThisRound=join([tifNamesThisRound;repmat({num2str(thisRound)},1,numChannels)],'',1);
        else
            idxScanField=ismember(formatSpec,namedScanFields);
            if sum(idxScanField)==0
                tifNamesThisRound=join([tifNamesThisRound;repmat({formatSpec},1,numChannels)],'',1);
            elseif sum(idxScanField)==1
                values=Scan.Rounds(thisRound).(formatSpec);
                
                if ~iscell(values)
                    if  isnumeric(values) && isvector(values) && length(values)==numChannels
                        values=strsplit(num2str(values)); % convert to cell array
                    else
                        error('makesStitches: tifNameFormat specifier %s failed in Scan round %i. %s should correspond to a cell array vector of same length as numChannels (%i)',formatSpec,thisRound,formatSpec,numChannels)
                    end
                elseif ~(isvector(values) && length(values)==numChannels)
                    error('makesStitches: tifNameFormat specifier %s failed in Scan round %i. %s should correspond to a cell array vector of same length as numChannels (%i)',formatSpec,thisRound,formatSpec,numChannels)
                end
                tifNamesThisRound=join([tifNamesThisRound;values],'',1);
            end
        end
    end
    % if any name ends with an underscore, then remove it
    for iChannel=1:numChannels
        tifName=tifNamesThisRound{iChannel};
        if endsWith(tifName,'_')
            tifNamesThisRound{iChannel}=tifName(1:end-1);
        end 
    end
    sTiffNames(thisRound).tifNames=tifNamesThisRound;
end

% check that the names are unique
allTiffNameVect={};
for thisRound=1:numRounds
   allTiffNameVect=[allTiffNameVect,sTiffNames(thisRound).tifNames];
end
numUniqueTifNames=length(unique(allTiffNameVect));
if numUniqueTifNames<length(allTiffNameVect)
    fprintf('based on tifNameFormat, these would be output tif file names. They are NOT unique. This is not allowed. Change tifNameFormat input to make them unique or do not input tifNameFormat to use defaults\n');
    fprintf('%s\n',allTiffNameVect{:})
    error('makesStitches: user-input tifNameFormat does not result in unique tif filenames')    
else
    fprintf('based on tifNameFormat, these will be the output tif file names:\n');
    fprintf('%s\n',allTiffNameVect{:})
end