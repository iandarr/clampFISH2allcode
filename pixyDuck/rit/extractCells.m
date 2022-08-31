function Tcell=extractCells(varargin)
% extracts cell-level data from the dataXXX.mat files
%
% Tcell=extractCells()
%   will look at all dataNNN.mat files in the current directory and extract the
%   cell (object) data assocated with areaCell, bbox, areaNucleus, isGood.
%   For example:
%
%   Tcell=
%
%     arrayNum    objNum    imageFileSuffix    areaCell    areaNucleus    isGood    bbox1    bbox2    bbox3    bbox4
%     ________    ______    _______________    ________    ___________    ______    _____    _____    _____    _____
%        1          1              1            37539          8160         1        205       59      225      381
%        1          2              1            37408          7541         1        218      475      203      457
%        2          1              2            43462          7445         1        710      739      301      264
%        2          2              2            12964          7188         1        631       10      182      119
%        3          1              4            32826         11705         1          3      446      324      157
%        3          2              4            36878          5680         1        437      505      402      437
%        3          3              4            39732          7720         1        769       10      240      355
%
%
%
% Tcell=extractCells(imageObjectTools)
%   will return one row of data for the one cell associated with object provided in imageObjectTools, instead of all objects in the current directory
%   you can get imageObjectTools by calling rajlabimagetools from within a directory with the data files:
%   Example:
%       tools=improc2.launchImageObjectTools;
%       tools.navigator.tryToGoToArray(2)
%       tools.navigator.tryToGoToObj(3)
%       Tcell=extractCells(tools)
%
% Tcell=extractCells(TextraMeasurements)  OR  Tcell=extractCells(imageObjectTools,TextraMeasurements)
%
%   will also add columns defined by TextraMeasurements. Example:
%
%     % TextraMeasurements =
%     %
%     %   5×4 table
%     %
%     %          columnName             area         valueToCompute    channel
%     %     ____________________    _____________    ______________    _______
%     %
%     %     'avgNuclearGFP'         'nuclear'          'average'        'gfp'
%     %     'avgCytoplasmGFP'       'cytoplasmic'      'average'        'gfp'
%     %     'avgTotalGFP'           'cell'             'average'        'gfp'
%     %     'stdTotalGFP'           'cell'             'std'            'gfp'
%     %     'avgExtracellularCY'    'nonCell'          'average'        'cy'
%     %
%     %
%     % columnName:
%     %   any character vector or string you choose, but must be a unique column name
%     %
%     % area:
%     %   'cell' is the full segmented area
%     %   'nonCell' is the area in the bounding box but outside the cell
%     %   'nuclear' is the automatic DAPI nuclear mask
%     %   'cytoplasmic' is 'cell' witth the 'nuclear' area removed
%     %
%     % valueToCompute:
%     %   'average' or 'avg' or 'mean' is the average (of the max merge)
%     %   'std' or 'standarddeviation' is the standard deviation (of the max merge)
%     %   'threshold' is not a calculation based on pixel values, but can be used for grabbing the threshold of a given channel
%     %
%     % channel:
%     %   choose 'gfp', 'tmr', 'cy', .... any channel you want
%
%
% Tcell=extractCells(...,'returnMask') will also output the binary mask of
% each object with the Tcell column variable name of cellMaskCropped
%



% check if imageObjectTools is given and set findObjectMode appropriately
toolsStruct=[];
TextraMeasurements=[];

myNargin=nargin;

% if 'returnMask' is an input, pull this out of varargin and update myNargin=myNargin-1
indexReturnMaskInput=strcmpi(varargin,'returnMasks');
if sum(indexReturnMaskInput)>=1
    returnMasks=true;
    varargin=varargin(~indexReturnMaskInput);
    myNargin=myNargin-sum(indexReturnMaskInput);
else
    returnMasks=false;
end

% if 'returnThreshold' is an input, pull this out of varargin and update myNargin=myNargin-2
isReturnThreshold=strcmpi(varargin,'returnThreshold');
if sum(isReturnThreshold)>1
    error('returnThreshold appears more than once in the input')
elseif sum(isReturnThreshold)==1
    returnThreshold=true;
    
    % look at next input
    indexReturnThreshold=find(isReturnThreshold);
    if length(varargin)<indexReturnThreshold+1
        error("the argument after 'returnThreshold' must be like {'cy', 'tmr'} or 'all', and currently there is no argument after 'returnThreshold'")
    end
    returnThresholdChannelsInput=varargin{indexReturnThreshold+1};
    if ischar(returnThresholdChannelsInput)
        returnThresholdChannelsInput={returnThresholdChannelsInput};
    elseif ~iscell(returnThresholdChannelsInput)
        error("argument after returnThreshold should be like {'cy','tmr'} or {'all'}")
    else % iscell
    end
    
    
    isReturnThresholdOrTheChannelsArgument=isReturnThreshold;
    isReturnThresholdOrTheChannelsArgument(indexReturnThreshold+1)=true;
    
    varargin=varargin(~isReturnThresholdOrTheChannelsArgument);
    myNargin=myNargin-2;
else
    returnThreshold=false;
end



if myNargin==2
    toolsStruct=varargin{1};
    TextraMeasurements=varargin{2};
    findObjectMode='providedAsInput';
    outputExtraMeasurements=true;
    assert(isstruct(toolsStruct))
    assert(istable(TextraMeasurements))
    
elseif myNargin==1
    if isstruct(varargin{1})
        toolsStruct=varargin{1};
        findObjectMode='providedAsInput';
        outputExtraMeasurements=false;
    elseif istable(varargin{1})
        TextraMeasurements=varargin{1};
        findObjectMode='findAllObjInCurrentDirectory';
        outputExtraMeasurements=true;
    else
        error('inputs should be tools structure or table')
    end
    
elseif myNargin==0
    findObjectMode='findAllObjInCurrentDirectory';
    outputExtraMeasurements=false;
else
    error('expecting zero, one, or two input arguments')
end


% check TextraMeasuremnents
if outputExtraMeasurements
    if ~isequal(TextraMeasurements.Properties.VariableNames,{'columnName','area','valueToCompute','channel'})
        error("expecting TextraMeasurments to have VariableNames 'columnName','area','valueToCompute','channel'")
    end
    
    extraMeasurementsColumnNames=[TextraMeasurements.columnName]';
    if ~all(ismember(lower(TextraMeasurements.area),lower({'nuclear','cytoplasmic','cell','nonCell'})))
        % total and cell are the same thing
        error(" Textrameasurements.area should be either 'nuclear','cytoplasmic','cell', or 'nonCell'")
    end
    
end

%% Extract from either just look at 1 object provided or look at all in this folder
% this calls the nested function extractCellFromToolsStruct
Tcell=table();
iCell=0;
switch findObjectMode
    case 'providedAsInput' % objectHandle is already a variable
        iCell=iCell+1;
        Tcell=extractCellFromToolsStruct;
        
    case 'findAllObjInCurrentDirectory' % get all objects in this directory
        toolsStruct = improc2.launchImageObjectTools;
        toolsStruct.iterator.goToFirstObject();
        while toolsStruct.iterator.continueIteration
            iCell=iCell+1;
            
            TcellThisObj=extractCellFromToolsStruct;
            if isempty(Tcell)
                Tcell=TcellThisObj;
            else
                Tcell=[Tcell;TcellThisObj];
            end
            
            toolsStruct.iterator.goToNextObject()
            
        end
    otherwise
        error("something messed up with inputs")
end


    function TcellThisObj=extractCellFromToolsStruct
        TcellThisObj=table();
        
        arrayNum=toolsStruct.navigator.currentArrayNum;
        objNum=toolsStruct.navigator.currentObjNum;
        
        % get file suffix (grab from the DAPI channel; if not dapi then the
        % first channel
        channelNames=toolsStruct.objectHandle.channelNames;
        if ismember('dapi',channelNames)
            imageFileName=toolsStruct.objectHandle.getImageFileName('dapi');
        else
            imageFileName=toolsStruct.objectHandle.getImageFileName(channelNames{1});
        end
        imageFileSuffix=str2num(char(join(regexp(imageFileName,'\d\d\d','match'),'')));
        
        % cell mask
        cellMaskCropped=toolsStruct.objectHandle.getCroppedMask;
        areaCell=sum(sum(cellMaskCropped));
        % nuclear mask
        nuclearMask=toolsStruct.objectHandle.getData('dapi').mask;
        %       NOTE: if it fails at the above line with error:
        %             Error using improc2.dataNodes.HandleToGraphBasedImageObject/findDataNode (line 261)
        %             Could not locate data of type improc2.interfaces.NodeData starting from node dapi.
        %             Error in improc2.dataNodes.HandleToGraphBasedImageObject/getData (line 70)
        %                   node = p.findDataNode(nodeLabel);
        %       Then you need to process image objects with improc2.processImageObjects before this
        
        areaNucleus=sum(sum(nuclearMask));
        
        bbox=toolsStruct.objectHandle.getBoundingBox;
        isGood=double(toolsStruct.annotations.getValue('isGood'));
        
        % store standard columns in TcellThisObj
        TcellThisObj.arrayNum(1)=arrayNum;
        TcellThisObj.objNum(1)=objNum;
        TcellThisObj.imageFileSuffix(1)=imageFileSuffix;
        TcellThisObj.areaCell(1)=areaCell;
        TcellThisObj.areaNucleus(1)=areaNucleus;
        TcellThisObj.isGood(1)=isGood;
        TcellThisObj.bbox1(1,:)=bbox(1);
        TcellThisObj.bbox2(1,:)=bbox(2);
        TcellThisObj.bbox3(1,:)=bbox(3);
        TcellThisObj.bbox4(1,:)=bbox(4);
        
        if returnMasks
            TcellThisObj.cellMaskCropped(1)={cellMaskCropped};
        end
        
        if returnThreshold
            objHandle=toolsStruct.objectHandle;
            labelsOfNodesWithData=objHandle.getLabelsOfNodesWithData;
            channelsWithSpots=replace(labelsOfNodesWithData(endsWith(labelsOfNodesWithData,':Spots')),':Spots','');
            
            if isequal(returnThresholdChannelsInput,{'all'})
                returnThresholdChannels=channelsWithSpots;
            else
                returnThresholdChannels=returnThresholdChannelsInput;
            end
            
            numberOfReturnThresholdChannels=length(returnThresholdChannels);
            
            for channelIndex=1:numberOfReturnThresholdChannels
                channelName=returnThresholdChannels{channelIndex};
                
                %% get :Spots node data (regionalMax data)
                thresholdInObj=objHandle.getData([channelName,':Spots']).threshold; % structure of this object and this channel's data, class aTrousRegionalMaxProcessedData
                
                TcellThisObj.([channelName,'Threshold'])(1)=thresholdInObj;
                %regionalMaxValues=spotDataStruct.regionalMaxValues;
                %regionalMaxIndices=spotDataStruct.regionalMaxIndices;
                %thresholdInObj=spotDataStruct.threshold;
            end
            
        end
        % add user-requested custom columns if they exist
        if outputExtraMeasurements
            
            for ii=1:height(TextraMeasurements)
                
                % get the right channel data (zMerge)
                % check that the channels requested are indeed in this
                % object, get that channel data
                channel=TextraMeasurements.channel{ii};
                if ismember(channel,channelNames)
                    %channelImg=toolsStruct.objectHandle.getData(channel).zMerge;
                    [zStack,~]=improc2.getCroppedImgAndMaskFromObj(toolsStruct.objectHandle,channel);
                    zMerge=max(zStack,[],3);
                else % that channel is not here. output nan's
                    warning(sprintf('the requested channel %s is not present, will output nan\n',channel));
                    zMerge=nan(size(cellMaskCropped,1),size(cellMaskCropped,2));
                end
                
                % get the right mask
                switch lower(TextraMeasurements.area{ii})
                    case 'nuclear'
                        mask=nuclearMask;
                    case 'cell'
                        mask=cellMaskCropped;
                    case 'noncell'
                        mask=~cellMaskCropped;
                    case 'cytoplasmic'
                        if isempty(nuclearMask)
                            mask=cellMaskCropped;
                        else
                            mask=cellMaskCropped&~nuclearMask;
                        end
                        
                    otherwise
                        error("the input area needing a custom value must be 'nuclear','cell','nonCell',or 'cytoplasmic'")
                end
                
                % compute value
                switch lower(TextraMeasurements.valueToCompute{ii})
                    case {'avg','average','mean'}
                        extraMeasurementVal=mean(zMerge(mask));
                    case {'std','standarddeviation'}
                        extraMeasurementVal=std(double(zMerge(mask)));
                    otherwise
                        error("valueToCompute is not acceptable")
                end
                
                
                % add computed value to table
                TcellThisObj.(TextraMeasurements.columnName{ii})(1)=extraMeasurementVal;
            end
        end
    end



end