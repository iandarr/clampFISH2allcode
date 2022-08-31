function excludeSlices(varargin)
% excludeSlices(numSlicesToExclude)
%   excludes the bottom numSlicesToExclude slices for all objects in the current directory
%
% excludeSlices(imageObjectBrowsingTools,numSlicesToExclude)
%   excludes the bottom numSlices slices, only for the object provided in imageObjectBrowsingTools
%   instead of all objects in the current directory
%   you can get this by calling rajlabimagetools from within a directory with the data files:
%   tools = improc2.launchImageObjectBrowsingTools();
%   tools.navigator.tryToGoToArray(1)
%   tools.navigator.tryToGoToObj(3)
%   excludeSlices(2, tools) will set the bottom 2 slices in array 1, object 3 to be excluded

if nargin==2
    if isstruct(varargin{1})
        if isequal(fieldnames(varargin{1}),{'navigator';'annotations';'objectHandle';'refresh'})
            findObjectMode='providedAsInput';
            tools=varargin{1};
            numSlicesToExclude=varargin{2};
        else
            error('if struct is provided as first input, must be tools from tools = improc2.launchImageObjectBrowsingTools()')
        end
    else
        error('first input must be tools struct from tools = improc2.launchImageObjectBrowsingTools()')
    end
elseif nargin==1
    findObjectMode='findAllObjInCurrentDirectory';
    numSlicesToExclude=varargin{1};
else
    error('could not interpret input')
end

assert(isnumeric(numSlicesToExclude))


switch findObjectMode
    case 'providedAsInput'
        % objectHandle is already a variable
        excludeBottomSlicesFunct(tools,numSlicesToExclude)
        tools.navigator.tryToGoToArray(1) %only done to save
    case 'findAllObjInCurrentDirectory'
        % Start rajlabimagetools in this directory
        tools = improc2.launchImageObjectBrowsingTools();
        
        for iArray=1:tools.navigator.numberOfArrays
            tools.navigator.tryToGoToArray(iArray)
            for iObj=1: tools.navigator.numberOfObjectsInCurrentArray
                tools.navigator.tryToGoToObj(iObj)
                
                excludeBottomSlicesFunct(tools,numSlicesToExclude)
            end
        end
        tools.navigator.tryToGoToNextObj %only done to save last object
    otherwise
        error("something messed up with inputs")
end

% nested function
    function excludeBottomSlicesFunct(tools,numBottomSlicesToExclude)
        objectHandle=tools.objectHandle;
        %excludeSlicesFromObject(objectHandle,numBottomSlicesToExclude)
        [rnaChannels, rnaProcessorClassName] = improc2.thresholdGUI.findRNAChannels(objectHandle);
        rnaChannelSwitch = dentist.utils.ChannelSwitchCoordinator(rnaChannels);
        rnaProcessorDataHolder = improc2.utils.ProcessorDataHolder(...
            objectHandle, rnaChannelSwitch, rnaProcessorClassName);
        sliceExcluder=improc2.utils.SliceExcluderForRegionalMaxProcData(rnaProcessorDataHolder);
        for iChannel=1:length(rnaChannels)
            channelName=rnaChannels{iChannel};
            rnaChannelSwitch.setChannelName(channelName);
            sliceExcluder.clearExclusionsAndExcludeSlicesUpTo(numBottomSlicesToExclude) % numBottomSlicesToExclud=3 here means spots from slice 1,2, and 3 are not counted, all other slices are counted
        end
        
    end






end