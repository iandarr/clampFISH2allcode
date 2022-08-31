function fitSpots(channels,varargin)
%  fitSpots(channels)
% 
% fitSpots(channels)
%   calls TwoStageSpotFitProcessedData for the input channels for all objects in current directory

if ischar(channels)
    channels={channels};
end
assert(size(channels,1)==1)
% if nargin>2
%     if strcmp(varargin{2},'OverwriteFittedData')
%         deleteFittedNodesOfDataFilesForChannels(channels)
%     else
%         error("third argument, if provided, must should be 'OverwriteFittedData'")
%     end
% end

dataAdder = improc2.processing.DataAdder();
unprocessedFittedData = improc2.nodeProcs.TwoStageSpotFitProcessedData();

for iChannel=1:length(channels)
    channelToFit=channels{iChannel};
    dataAdder.addDataToObject(unprocessedFittedData, channelToFit, [channelToFit,':Fitted'])
end

dataAdder.repeatForAllObjectsAndQuit();
improc2.processing.updateAll 

end



% 
% function deleteFittedNodesOfDataFilesForChannels(channelsForWhichToDeleteFittedData)
% channelsForWhichToDeleteFittedData=join([channelsForWhichToDeleteFittedData',repmat({':Fitted'},length(channelsForWhichToDeleteFittedData),1)],'',2)';
% 
% dataFilenamesStruct=dir('data*.mat');
% if isempty(dataFilenamesStruct)
%     error('no rajlabimagetools dataXXX.mat files found in directory %s',pwd)
% else
%     dataFilenames=[{dataFilenamesStruct.name}];
% 
%     for  iFile=1:length(dataFilenames)
%         dataFilename=dataFilenames{iFile};
%         load(dataFilename,'objects')
%         if isempty(objects)
%             % moving on to next dataXXX.mat file
%         else
% 
%             numObjects=length(objects);
% 
%             for iObj=1:numObjects
%                 thisObject=objects(iObj);
%                 allNodeNames=arrayfun(@(x) {x{1}.label}, thisObject.graph.nodes);
%                 idxNodesToDelete=ismember(allNodeNames,channelsForWhichToDeleteFittedData);
%                 if ~any(idxNodesToDelete)
%                     fprintf('  For object number %i no nodes in %s have any fitted data for these channels\n',iObj,dataFilename)
%                 else
%                     % remove these nodes
%                     thisObject.graph.nodes=thisObject.graph.nodes(~idxNodesToDelete);
%                     %thisObject.graph.labels=thisObject.graph.labels(~idxNodesToDelete);
%                     thisObject.graph.childVsParentConnectivity=thisObject.graph.childVsParentConnectivity(~idxNodesToDelete,~idxNodesToDelete);
%                     objects(iObj)=thisObject;
%                 end
%             end
%             save(dataFilename,'objects')
%         end
%     end
% end
% 
% end