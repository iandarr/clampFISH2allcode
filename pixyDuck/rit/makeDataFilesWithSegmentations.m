function varargout=makeDataFilesWithSegmentations(destinationDir,fileNumList,maskImgList,varargin)
% makeDataFilesWithSegmentations(destinationDir,filenumList,maskImgList)
%
% useful for:
%      - creating segmentation data files from automatically-generated
%      segmentations (Ie. from cellpose)
%      - copying segmentations from one frame of reference (Ie. 60X) into another
%        destination frame of reference (Ie. 10X), as long as the destination
%        fields of view actually cover that area.
%
% fileNumList is a numeric column vector
% 
% maskImgList is a cell array (column) where each element is a binary mask image (true or 1 = cell, false or 0 =not cell). Doesn't have to be logical class.
%
% destinationDir is a directory name to output the dataXXX.mat files
%% Mock inputs before I make this a function
% maskImg=false(2044,2048);
% maskImg(10:100,500:1500)=true;
% imshow(maskImg)
% imagenumber='001';
% destinationDir='/Volumes/IAND_04/data/E157_validation/VsSmFISH/20X/s13_PreStrip_PresmFISH/w1256/Tiff_s13_c2w2_naive_EGFR';
%
%
p=inputParser;
p.addParameter('ignoreSmallNonconnectedSegmentations',false,@islogical); % if the segmentation mask has a non-connecting region (Ie. an 'island'), then ignore this
p.parse(varargin{:});
ignoreSmallNonconnectedSegmentations=p.Results.ignoreSmallNonconnectedSegmentations;

assert(length(fileNumList)==length(maskImgList))

indOfSegmentationsToOutput=find(~isnan(fileNumList));
if isempty(indOfSegmentationsToOutput)
    warning('all %i rows of fileNumList are nan. This can occur if all segmentations are out of bounds of destination fields of view',length(fileNumList))
end
numSegmentationsToOutput=length(indOfSegmentationsToOutput);

if ~isdir(destinationDir)
   error('could not find destinationDir=%s',destinationDir) 
end

%%  first delete data00X.mat files in the destination folder
dataFilesInDestinationDir=listFilesWithExtension('.mat','directoryPath',destinationDir);
dataFilesInDestinationDir=dataFilesInDestinationDir(startsWith(dataFilesInDestinationDir,'data'));

for iFileToDelete=1:length(dataFilesInDestinationDir)
    delete(fullfile(destinationDir,filesep,dataFilesInDestinationDir{iFileToDelete}))
end
%delete(sprintf('%s%sdata%03d.mat',destinationDir,filesep,fileNum)

%% output blank data files with objects=[] in destination for all images
[~,fileNums,~]=getImageFiles(destinationDir);
for fileNum=fileNums
    objects=[];
    %fileNumStr=sprintf('%03d',fileNum);
    save(sprintf('%s%sdata%03d.mat',destinationDir,filesep,fileNum),'objects');
end

%% Add segmentations to blank data00X.mat files. Loop through objects in list
for iSegmentationOut=1:numSegmentationsToOutput
    ind=indOfSegmentationsToOutput(iSegmentationOut);
    fileNum=fileNumList(ind);
    fileNumStr=sprintf('%03d',fileNum);
    
    maskImg=maskImgList{ind};
    
    % check the mask for multiple regions
    if ignoreSmallNonconnectedSegmentations
        CC=bwconncomp(maskImg);
        if CC.NumObjects>1
            % then we have 2 or more objects
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggestNumPixels,idxBiggest] = max(numPixels);
            maskImg=false(size(maskImg));
            maskImg(CC.PixelIdxList{idxBiggest})=true; % only keep the biggest object's mask
            
%             [smallestNumPixels,idxSmallest] = min(numPixels);
%             maskImgSmallest=false(size(maskImg));
%             maskImgSmallest(CC.PixelIdxList{idxSmallest})=true;
            
        end
            
    end
    % this is called in improc2.segmentGUI.SegmentGUI under function segmentObject_Callback
    newObj=improc2.buildImageObject(maskImg,fileNumStr,destinationDir);
    
    % open the data file and re-save with new newObj added
    loadData=load(sprintf('%s%sdata%03d.mat',destinationDir,filesep,fileNum));
    objects=loadData.objects;
    objects=[objects,newObj];
    
    save(sprintf('%s%sdata%03d.mat',destinationDir,filesep,fileNum),'objects');
    
end

if nargout==2
    % varargout{1} is dArrayNum
    % varargout{2} is dObjNum
    dArrayNum=fileNumList; % these two variables are equivalent only if there are is a data00X.mat file for every single field of view, even if there are no objects segmented in some of them.
    
    dObjNum=nan(length(fileNumList),1);
    uniqueFileNumListRowvect=unique(fileNumList)';
    uniqueFileNumListRowvect=uniqueFileNumListRowvect(~isnan(uniqueFileNumListRowvect));
    for thisFileNum=uniqueFileNumListRowvect
        indsOfThisFileNum=find(fileNumList==thisFileNum);
        for iObjInThisArray=1:length(indsOfThisFileNum)
            ind=indsOfThisFileNum(iObjInThisArray);
            dObjNum(ind)=iObjInThisArray;
        end
    end
    
    varargout=cell(1,2);
    varargout{1}=dArrayNum;
    varargout{2}=dObjNum;
end
end