function setThreshold(channels,thresholds,varargin)
% setThreshold(channels,thresholds)
%   sets thresholds for all objects in current directory
%
% setThreshold(channels,thresholds,imageObjectTools)
%   sets threshold for only the object provided by imageObjectTools
%   instead of all objects in the current directory
%   you can get this by calling rajlabimagetools from within a directory with the data files:
% 
%
%   Example 1 (set threshold for all objects in current directory, for tmr and cy channels)
%       setThreshold({'tmr','cy'},[200 250]) % tmr threshold is set ot 200, cy threshold is set to 250
%
%   Exampld 2 (set threshold for all objects in current directory, for tmr channel only)
%
%       setThreshold('tmr',200) % or, equivalently:
%       setThreshold({'tmr'},200)
% 
%   Example 3 (set threshold for only Array 2, Object 3, channels tmr and cy):       
%       tools=improc2.launchImageObjectTools;
%       tools.navigator.tryToGoToArray(2)
%       tools.navigator.tryToGoToObj(3)
%       setThreshold({'tmr','cy'},[200 250], tools)
%


if nargin==3
    if isstruct(varargin{1})
        tools=varargin{1};
        findObjectMode='providedAsInput';
    else
        error('inputs should be tools structure or table')
    end
else
    findObjectMode='findAllObjInCurrentDirectory';
end

if ischar(channels)
    channels={channels};
end

if ~isnumeric(thresholds)
   error("thresholds must be a matrix of a size that matches the cell array channels")
end

if ~isequal(size(channels),size(thresholds))
    error("channels and thresholds must be the same size. Put channels into a cell array of character vectors and thresholds into a matrix")
end



switch findObjectMode
    case 'providedAsInput'
        % objectHandle is already a variable
        setThresholdFunct
        
        tools.navigator.tryToGoToArray(1) %only done to save
    case 'findAllObjInCurrentDirectory'
        % Start rajlabimagetools in this directory
        tools = improc2.launchImageObjectBrowsingTools();
        
        for iArray=1:tools.navigator.numberOfArrays
            tools.navigator.tryToGoToArray(iArray)
            for iObj=1: tools.navigator.numberOfObjectsInCurrentArray
                tools.navigator.tryToGoToObj(iObj)
                
                setThresholdFunct
            end
        end
        tools.navigator.tryToGoToNextObj %only done to save last object
    otherwise
        error("something messed up with inputs")
end


    function setThresholdFunct % nested function
        for iChannel=1:length(channels)
            channel=channels{iChannel};
            threshold=thresholds(iChannel);
            
            dataToBeModified=tools.objectHandle.getData([channel,':Spots']);
            dataToBeModified.threshold=threshold;
            %dataToBeModified.needsUpdate=true;
            tools.objectHandle.setData(dataToBeModified, [channel,':Spots'])
            dataToBeModified=tools.objectHandle.getData([channel,':Spots']);
            
            %tools.navigator.saveIfNeedsSave
        end
        
    end



end