function excludeSlicesFromObject(objectHandle,numBottomSlicesToExclude)
% will call clearExclusionsAndExcludeSlicesUpTo
% a method of improc2.utils.SliceExcluderForRegionalMaxProcData
% objectHandle can be gotten by:
%   browsingTools = improc2.launchImageObjectBrowsingTools();
%   objectHandle = browsingTools.objectHandle;
%
% you can use the navigator likes this:
%   browsingTools.navigator.goToArray(1)
% 


[rnaChannels, rnaProcessorClassName] = improc2.thresholdGUI.findRNAChannels(objectHandle);

rnaChannelSwitch = dentist.utils.ChannelSwitchCoordinator(rnaChannels);
    
rnaProcessorDataHolder = improc2.utils.ProcessorDataHolder(...
        objectHandle, rnaChannelSwitch, rnaProcessorClassName);

%rnaProcessorDataHolder.processorData; %this field exists

sliceExcluder=improc2.utils.SliceExcluderForRegionalMaxProcData(rnaProcessorDataHolder);
sliceExcluder.clearExclusionsAndExcludeSlicesUpTo(numBottomSlicesToExclude) % numBottomSlicesToExclud=3 here means spots from slice 1,2, and 3 are not counted, all other slices are counted

end