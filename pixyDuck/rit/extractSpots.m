function TspotsOut=extractSpots(channels,varargin)
% TspotsOut=extractSpots(channelList)
%   will return only spot intensities above the object's threshold value
%   like this:
%      arrayNum         objNum       channel     X      Y     Z    intensities    isGood    threshold
%     ___________    ____________    _______    ___    ___    _    ___________    ______    _________
%          1              2           'tmr'     116    220    1      37.006       true         45
%          1              2           'tmr'     130    101    5      45.091       true         45
%          1              2           'tmr'     127    160    3      45.121       true         45
%          1              3           'tmr'     136    197    4      35.123       true         100.45
%          1              3           'tmr'      41    370    1      200.13       true         100.45
%
%   and if those spots have an associated :Fitted node (from Guassian
%   fitting), this will also include addition columns with fitted data:
%
%
% TspotsOut=extractSpots(channelList,'allSpotIntensities')
%   will return all spot intensities of all objects in the current directory
%   even if those intensities are below the object's threshold value
%
% TspotsOut=extractSpots(channelList,'intensitiesAboveThreshold')
%   will return only spots above the threshold value currently stored with the object
%
% TspotsOut=extractSpots(channelList,'intensitiesAboveGivenCutoff',cutoffValues)
%   returns only intensities equal to or above a user-supplied
%   cutoffValue(s), ignoring the threshold stored in the object. Must
%   provide one cutoff value for each channel in channelList. If
%   channelList is 'all', then only length(thresholdValues) should be 1,
%   and all channels will use this single cutoffValue.
%
% TspotsOut=extractSpots(channelList,imageObjectTools,...)
%   will return spots only for the object provided in imageObjectTools,
%   instead of all objects in the current directory
%   you can get this by calling rajlabimagetools from within a directory with the data files:
%   Example:       
%       tools=improc2.launchImageObjectTools;
%       tools.navigator.tryToGoToArray(2)
%       tools.navigator.tryToGoToObj(3)
%       TspotsOut=extractSpots({'tmr', 'cy'}, tools)
%
% TspotsOut=extractSpots(channelList,...,'getFittedData',true)
%   in addition to normal spot XY and intensity data (from
%   RegionalMaxValues), will also look for a :Fitted node from Gaussian
%   fitting and, if that data is present, return it in 4 additional
%   columns, like this:
%
%
%
% channels is a cell array of channels you want to extract:
%       channels={'tmr'};
%       channels={'trans','dapi','gfp','tmr','alexa','cy','nir'};
%
%   OR it can be just a character vector of one of these channels:
%   channels='tmr'
%
% 	OR it channels can be 'all' which returns all the object's spot data
%

if ischar(channels)
    channels=cellstr(channels);
end

if ~iscell(channels)
    error('channelList should be a cell array of character vectors, or a character vector for a single channel')
end


%cutoffValueOrValues=nan(1,length(channels));


% figure out if toolsStruct has been provided. take it from varargin and
% then make vararginWithoutToolsStruct
if nargin>=2
    isStructVarargin=false(1,length(varargin));
    for iVarargin=1:length(varargin)
        isStructVarargin(iVarargin)=isstruct(varargin{iVarargin});
    end

    iVararginIsToolsStruct=find(isStructVarargin);
    if ~(length(iVararginIsToolsStruct)<=1)
        error('found more than 1 struct input')
    end
    
    if ~isempty(iVararginIsToolsStruct)
        toolsStruct=varargin{iVararginIsToolsStruct};
        findObjectMode='providedAsInput';
        vararginRemaining=varargin(~isStructVarargin);
    else
        findObjectMode='findAllObjInCurrentDirectory';
        vararginRemaining=varargin;
    end
else
    vararginRemaining=varargin;
    findObjectMode='findAllObjInCurrentDirectory';
end


% intensityMode (what spots to extract)
if isempty(vararginRemaining)
    intensityMode='intensitiesAboveThreshold';
else
    if strcmpi(vararginRemaining{1},'allSpotIntensities')
        intensityMode='allSpotIntensities';
        if length(vararginRemaining)>1,
            vararginRemaining=vararginRemaining(2:end);
        else
            vararginRemaining={};
        end
    elseif strcmpi(vararginRemaining{1},'intensitiesAboveThreshold')
        intensityMode='intensitiesAboveThreshold';
        if length(vararginRemaining)>1,
            vararginRemaining=vararginRemaining(2:end);
        else
            vararginRemaining={};
        end
    elseif strcmpi(vararginRemaining{1},'intensitiesAboveGivenCutoff')
        intensityMode='intensitiesAboveGivenCutoff';
        if length(vararginRemaining)<2
            error("if argument is intensitiesAboveGivenCutoff then the next argument (cutoffValue) must be a number")
        end
        if ~isnumeric(vararginRemaining{2})
            error("third argument (cutoffValue) must be a number or numbers")
        else
            cutoffValueOrValues=vararginRemaining{2};
            
            % check that user-provided cutoffValueOrValues is consistent
            % with the channels to extract:
                % if 3 channels to extract, should have 3 cutoff values
                % if 'all' channels, then only 1 cutoff value is allowed
           if isequal(channels,{'all'})
               if length(cutoffValueOrValues)~=length(channels)
                   error(sprintf("if channels input is 'all' with 'intensitiesAboveGivenCutoff' true, then input cutoffValues must be a length of 1"))
               end
           elseif length(cutoffValueOrValues)~=length(channels)
               error(sprintf("if you extract %i channels with 'intensitiesAboveGivenCutoff' true, then the input cutoffValues must be a 1 x %i numeric matrix", length(channels),length(channels)))

           end
        end
        
        if length(vararginRemaining)>2
            vararginRemaining=vararginRemaining(3:end);
        else
            vararginRemaining={};
        end
    else
        error("first non-imageObjectTools argument must be allSpotIntensities, intensitiesAboveThrehold, or intensitiesAboveGivenCutoff")
    end
end




%% read any remaining inputs (besides tools channelList, imageObjectTools, and intensityMode

p=inputParser;
p.addParameter('getFittedData',false,@islogical)
p.parse(vararginRemaining{:});
getFittedData=p.Results.getFittedData;

% check channels
if ~isequal(channels,{'all'}) % this is okay
    acceptableChannels={'trans','dapi','gfp','gfpa','gfpb','gfpc','gfpd','gfpe','tmr','tmra','tmrb','tmrc','tmrd','tmre','alexa','alexaa','alexab','alexac','alexad','alexae','cy','cya','cyb','cyc','cyd','cye','nir','gfpL','tmrL','alexaL','cyL'};
    if(~all(ismember(channels,    acceptableChannels)))
        error('at least one of the input channels not accepted. Consider editing the acceptableChannels variable within this function to reflect your channel names')
    end
    
end




%% Extract from either just look at 1 object provided or look at all in this folder
% this calls the nested function extractSpotsFromToolsStruct
TspotsOut=[];
switch findObjectMode
    case 'providedAsInput'
        % objectHandle is already a variable
        TspotsOut=extractSpotsFromToolsStruct;
        
    case 'findAllObjInCurrentDirectory'
        % Start rajlabimagetools in this directory
        toolsStruct = improc2.launchImageObjectTools;
        toolsStruct.iterator.goToFirstObject();
        while toolsStruct.iterator.continueIteration
            
            
            TspotsThisObj=extractSpotsFromToolsStruct;
            if isempty(TspotsOut)
                TspotsOut=TspotsThisObj;
            else
                TspotsOut=[TspotsOut;TspotsThisObj];
            end
            
            toolsStruct.iterator.goToNextObject()
        end
    otherwise
        error("something messed up with inputs")
end

%% nested function for extraction from an object
    function TspotsThisObj=extractSpotsFromToolsStruct
        
        
        objHandle=toolsStruct.objectHandle;
        labelsOfNodesWithData=objHandle.getLabelsOfNodesWithData;
        channelsWithSpots=replace(labelsOfNodesWithData(endsWith(labelsOfNodesWithData,':Spots')),':Spots','');
        
        if isequal(channels,{'all'})
            channelListToExtract=channelsWithSpots;
        else
            channelListToExtract=channels;
        end
        
        % channels with :Spots node may or may not have :Fitted node. Get a
        % list of them to extract later on
        channelsWithFitted=replace(labelsOfNodesWithData(endsWith(labelsOfNodesWithData,':Fitted')),':Fitted','');
        
        if ~all(ismember(channelListToExtract,channelsWithSpots))
            error("not all the channels you are trying to extract have spots, or they may not be a channel at all")
        end
        
        numberOfChannelsToExtract=length(channelListToExtract);
        % check over these
        
        for channelIndex=1:numberOfChannelsToExtract
            channelName=channelListToExtract{channelIndex};
            
            %% get :Spots node data (regionalMax data)
            spotDataStruct=objHandle.getData([channelName,':Spots']); % structure of this object and this channel's data, class aTrousRegionalMaxProcessedData
            
            regionalMaxValues=spotDataStruct.regionalMaxValues;
            regionalMaxIndices=spotDataStruct.regionalMaxIndices;
            thresholdInObj=spotDataStruct.threshold;
            
            switch intensityMode
                case 'intensitiesAboveThreshold'
                    spotInds=regionalMaxValues>thresholdInObj;
                    %[X, Y, Z] = spotDataStruct.getSpotCoordinates();
                    %intensities = regionalMaxValues(regionalMaxValues>thresholdInObj);
                    %for checking over
                    %[Xcheck,Ycheck,Zcheck]=ind2sub(spotDataStruct.imageSize,regionalMaxIndices(regionalMaxValues>thresholdInObj));
                    %assert(all([X==Xcheck;Y==Ycheck;Z==Zcheck],1))
                case 'intensitiesAboveGivenCutoff'
                    if isequal(channels,{'all'})
                        cutoffValue=cutoffValueOrValues;
                    else
                        cutoffValue=cutoffValueOrValues(channelIndex);
                    end
                    spotInds=regionalMaxValues>cutoffValue;
                    %[X,Y,Z]=ind2sub(spotDataStruct.imageSize,regionalMaxIndices(regionalMaxValues>cutoffValue));
                    %intensities = regionalMaxValues(regionalMaxValues>thresholdInObj);
                case 'allSpotIntensities'
                    spotInds=true(size(regionalMaxValues));
                otherwise
                    error("intensityMode should be 1 of 3 cases")
            end
            
            
            [Y,X,Z]=ind2sub(spotDataStruct.imageSize,regionalMaxIndices(spotInds)); %if there are excluded slices this may not work
            intensities=regionalMaxValues(spotInds);
            
            % check the 'above threshold case since we can get the
            % values by the method getSpotCoordinates()
            
            switch intensityMode
                case 'intensitiesAboveThreshold'
                    [Ycheck,Xcheck,Zcheck]=spotDataStruct.getSpotCoordinates();
                    assert(all([X==Xcheck;Y==Ycheck;Z==Zcheck],1))
                    Xcheck=[];Ycheck=[];Zcheck=[];
            end
            
            
            %% get object-level data
            numSpotsGrabbed=sum(spotInds);
            assert(numSpotsGrabbed==size(X,1));
            
            objNum=toolsStruct.navigator.currentObjNum();
            objectNumberColVect = repmat(objNum, [numSpotsGrabbed 1]);
            arrayNum=toolsStruct.navigator.currentArrayNum();
            arrayNumberColVect = repmat(arrayNum, [numSpotsGrabbed 1]);
            channelColVect = repmat({channelName}, [length(X) 1]);
            
            isGoodColVect = repmat(toolsStruct.annotations.getValue('isGood'), [numSpotsGrabbed 1]);
            
            thresholdColVect = repmat(thresholdInObj, [numSpotsGrabbed 1]);
            
            %% Add all :Spots data to a table
            TspotsThisObjThisChannel=[array2table(nan(numSpotsGrabbed,7),'VariableNames',{'arrayNum','objNum','channel','X','Y','Z','intensities'}),...
                table('Size',[numSpotsGrabbed 1],'VariableTypes',{'logical'},'VariableNames',{'isGood'}),...
                array2table(nan(numSpotsGrabbed,1),'VariableNames',{'threshold'})];
            
            
            TspotsThisObjThisChannel.objNum=objectNumberColVect;
            TspotsThisObjThisChannel.arrayNum=arrayNumberColVect;
            TspotsThisObjThisChannel.channel=channelColVect;
            TspotsThisObjThisChannel.X=X;
            TspotsThisObjThisChannel.Y=Y;
            TspotsThisObjThisChannel.Z=Z;
            TspotsThisObjThisChannel.intensities=intensities;
            TspotsThisObjThisChannel.isGood=isGoodColVect;
            TspotsThisObjThisChannel.threshold=thresholdColVect;
            
            %% get :Fitted node data (from gaussian fitting) if it exists
            if getFittedData
                if ~ismember(channelName,channelsWithFitted)
                    % then this channel does not have :Fitted data
                    % pad with nan.
                    TspotsThisObjThisChannel=[TspotsThisObjThisChannel,...
                        array2table(nan(height(TspotsThisObjThisChannel),4),'VariableNames',{'xFitted','yFitted','amplitudeFitted','sigmaFitted'})];
                else
                    % then this object and channel has :Fitted data, so
                    % get Fitted spots. Spot fitting is only done for spots
                    % with RegionalMax values above the object's stored threshold
                    % value, so this function may return more or less spots
                    % than are in above TspotsThisObjThisChannel table
                    
                    fittedSpotDataStruct=objHandle.getData([channelName,':Fitted']).getFittedSpots();
                    numFittedSpots=length(fittedSpotDataStruct);
                    
                    if numFittedSpots==0
                        TfittedSpotsThisObjThisChannel=array2table(nan(numSpotsGrabbed,4),...
                            'VariableNames',{'xFitted','yFitted','amplitudeFitted','sigmaFitted'});
                        
                        TspotsThisObjThisChannel=[TspotsThisObjThisChannel,TfittedSpotsThisObjThisChannel];
                        
                    else % then there are Fitted Spots - get these, then match these up with the :Spots data in TspotsThisObjThisChannel
                        TfittedSpotsThisObjThisChannel=array2table([...
                            [fittedSpotDataStruct.xCenter]',...
                            [fittedSpotDataStruct.yCenter]',...
                            [fittedSpotDataStruct.zPlane]',...
                            [fittedSpotDataStruct.amplitude]',...
                            [fittedSpotDataStruct.sigma]'],...
                            'VariableNames',{'xFitted','yFitted','zFitPlaneDeleteMe','amplitudeFitted','sigmaFitted'});
                        
                        
                        % check that the number of fitted spots (wich are only for above-threshold regionalMax values) are consistent with the
                        % current threshold in the object.
                        numAllObjSpotsAboveThreshold=sum(regionalMaxValues>thresholdInObj);
                        if numAllObjSpotsAboveThreshold~=numFittedSpots
                            warning('Array %i Obj %i, %s: there are %i spots above current threshold (%g) but %i fitted spots. may need to re-do fitting with this updated threshold value',arrayNum,objNum,channelName,numAllObjSpotsAboveThreshold,thresholdInObj,numFittedSpots)
                        end
                        
                        fprintf('Array %i Obj %i, %s: grabbing %i spots. There are %i spots above current threshold (%g) and %i fitted spots\n',arrayNum,objNum,channelName,numSpotsGrabbed,numAllObjSpotsAboveThreshold,thresholdInObj,numFittedSpots)
                        
                        % now add :Fitted data to the :Spots data. First, need
                        % to adjust the number of rows if there is a mismatch
                        % between number of fitted spots (numFittedSpots) and
                        % number of regionalMax values obtained (numSpotsGrabbed)
                        
                        %numGrabbedSpotsAboveThreshold=sum(intensities>=thresholdColVect);
                        
                        if numSpotsGrabbed>numFittedSpots
                            
                            %then pad :Fitted data with nan
                            numRowsToNanPad=numSpotsGrabbed-numFittedSpots;
                            colVariableNames=TfittedSpotsThisObjThisChannel.Properties.VariableNames;
                            numColsToPad=length(colVariableNames);
                            TfittedSpotsThisObjThisChannel=[array2table(nan(numRowsToNanPad,numColsToPad),'VariableNames',colVariableNames);TfittedSpotsThisObjThisChannel];
                        elseif numSpotsGrabbed<numFittedSpots
                            TfittedSpotsThisObjThisChannel=TfittedSpotsThisObjThisChannel(end-numSpotsGrabbed+1:end,:);
                        end
                        
                        % merge the two tables
                        TspotsThisObjThisChannel=[TspotsThisObjThisChannel,TfittedSpotsThisObjThisChannel];
                        
                        % check that Z values are identical - they should be. If they are,
                        % delete this column
                        if ~all(any([isnan(TspotsThisObjThisChannel.zFitPlaneDeleteMe),...
                                TspotsThisObjThisChannel.Z==TspotsThisObjThisChannel.zFitPlaneDeleteMe],2))
                            error("Array %i, Obj %i, %s: mismatch between :Spots and :Fitted data. Resolve this (Can occur if excludeSlices was different when spot fitting was performed). Or run with 'getFittedData',false",arrayNum,objNum,channelName)
                        else
                            TspotsThisObjThisChannel.zFitPlaneDeleteMe=[];
                        end
                        
                        % check that the X and Y values are close
                        
                        distErrorXY=(TspotsThisObjThisChannel.X - TspotsThisObjThisChannel.xFitted).^2 + (TspotsThisObjThisChannel.Y - TspotsThisObjThisChannel.yFitted).^2;
                        if ~all(any([distErrorXY<4,isnan(distErrorXY)],2))
                            warning('Array %i, Obj %i, %s: at least one spot from :Spots node (regionalMaxValues) and :Fitted node has significantly different XY coordinates. Can occur if excludeSlices is different from when spot fitting was performed. Or if allowable range for gaussian fit is large',arrayNum,objNum,channelName)
                        end
                        
                    end
                    
                end
            end % end "if getFittedData"
            
            %% Add to table of other objects
            if channelIndex==1
                TspotsThisObj=TspotsThisObjThisChannel;
            else
                TspotsThisObj=[TspotsThisObj;TspotsThisObjThisChannel];
            end
        end % end channel loop
        
    end % end nested function extractSpotsFromObjHandle

end %end main extractSpots function


