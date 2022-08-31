function [imgCropped,croppedMask,spotRowColZ]=getRajlabimagetoolsCroppedCellImage(imgDir,arrayNum,objNum,channel,spotThreshold,zplanesInput,getSpots)

% imshow(toolsThresh.browsingTools.objectHandle.getData('dapi').zMerge)

%toolsThresh=improc2.launchThresholdGUI;
%toolsThresh.rnaChannelSwitch.setChannelName(sourceChannel)
%toolsThresh.browsingTools.navigator.tryToGoToArray(sourceArrayNum)
%toolsThresh.browsingTools.navigator.tryToGoToObj(sourceObjNum)

origDir=pwd;

cd(imgDir)
tools=improc2.launchImageObjectTools;
cd(origDir);


tools.navigator.tryToGoToArray(arrayNum)
tools.navigator.tryToGoToObj(objNum)


% get mask
mask=tools.objectHandle.getMask;

% cropped mask
croppedMask=tools.objectHandle.getCroppedMask;


% get spots
if getSpots
% threshold
if iscell(spotThreshold)
    spotThreshold=spotThreshold{1};
end
if ischar(spotThreshold)&&strcmp(spotThreshold,'manual')
    TspotsOut=extractSpots(channel,tools,'intensitiesAboveThreshold');
elseif isnumeric(spotThreshold)
    TspotsOut=extractSpots(channel,tools,'intensitiesAboveGivenCutoff',spotThreshold);
else error('huh')
end
spotRowColZ=[TspotsOut.Y,TspotsOut.X,TspotsOut.Z];

else
    spotRowColZ=[];
end
%isSpot=regionalMaxima>=spotThreshold;
%spotIntensities=tools.objectHandle.getData(channel).regionalMaxValues(isSpot);

%spotIdx=tools.objectHandle.getData(channel).regionalMaxIndices(isSpot);

%imgSize=tools.objectHandle.getData(channel).imageSize;

%[spotRows,spotCols,spotsZ]=ind2sub(imgSize,spotIdx);
%spotRowColZ=[spotRows,spotCols,spotsZ];


%imgPath=fullfile(imgDir,[channel,sprintf('%03i',arrayNum),'.tif']);
imgFileName=tools.objectHandle.getImageFileName(channel);
imgPath=fullfile(imgDir,filesep,imgFileName);
            numZplanesInImage=length(imfinfo(imgPath));
            if strcmp(zplanesInput,'all')
                zplanesToGet=1:numZplanesInImage;
            else
                zplanesToGet=zplanesInput;
            end

            for zplane=zplanesToGet
                if zplane==zplanesToGet(1) % first plane
                    img=imread(imgPath,zplane);
                else
                    img=max(cat(3,img,imread(imgPath,zplane)),[],3);
                end
            end

% bounding box
bbox=tools.objectHandle.getBoundingBox;

imgCropped=img(bbox(2):bbox(2)+bbox(4),bbox(1):bbox(1)+bbox(3),:);

end