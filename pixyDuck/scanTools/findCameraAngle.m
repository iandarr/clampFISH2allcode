function [cameraAngle,cameraAngleList,numInliersList,statusFlagList]=findCameraAngle(ND2filenameORreader,stagecoor_XY,UmPerPixel,StageAndImageAxesAreMirrored,varargin)
% see findCameraOrientation.m for a more up-to-date function
% 
% cameraAngle=findCameraAngle(ND2filenameORreader,stagecoor_XY,UmPerPixel,StageYDimensionIsNotFlipped,estimateTformOptionalNameValueArguments)
%
% finds camera angle of an Nd2 file, looking at multiple adjacent pairs of
% images, taking into account their stage XY positions, and reaching a consensus (median) angle
%
% Required Inputs:
%   ND2filenameORreader: filepath to ND2 file or its reader
%
%   stagecoor_XY is two-column matrix of X (first column) and Y (second
%       column) stage positions with units of microns. Can get these with nd2positions()
%
%   UmPerPixel is microns per pixel of the image
%
%
%   StageYDimensionIsNotFlipped:
%       If you force the stagecoord +X dimension to point in the
%       same direction as the imagecoord +col dimension, then:
%            if  the stagecoord +Y dimension points in the +row dimension,
%                   then StageYDimensionIsNotFlipped=true (should be default)
%            if the stagecoord +Y dimension points in the -row dimension,
%                   then StageYDimensionIsNotFlipped=false
%       the same direction as the X coordinate, then is the Y d
%
%
% Optional Name-value inputs:
%   'tileIDpairs'
%       2-column vector of pairs of tileIDs to compare to one another
%       Example:
%           tileIDpairs=[1 2; 6 7; 19 20; 20 21; 21 22; 19 22]
%               note if the first row in the scan is 20 tiles long and then
%               snakes around, then tile 19 and 22 are next to eachother.
%
%   'numImagePairsToTryToFindAngleFor'
%       Example:
%           numImagePairsToTryToFindAngleFor=30;
%               will try to look at 30 (regularly spaced) image pairs and
%               output the median angle. Ignored if tileIDpairs is used
%   Name-Value pairs that will get passed down through
%       findCameraAngleBetweenTwoImages and then used by estimateTform. See
%       estimateTform for details
%
%
% OUTPUTS
% 
% statusFlag from estimateGeometricTransform2D:
%   0: No error.
%   1: matchedPoints1 and matchedPoints2 do not contain enough points.
%   2: Not enough inliers per threshold intrinsic to estimateGeometricTransform2D
% statusFlag from esimateTform:
%   3: not enough inliers per the threshold in estimateTform
% statusFlag from findCameraAngleBetweenTwoImages:
%   4: images are not close enough to eachtother
%% input parser

% for i=1:numel(varargin)
%     if  isempty(varargin{i})
%         keepIdx=true(1,numel(varargin)); keepIdx(i)=false; varargin=varargin(keepIdx);
%     end
% end

ip=inputParser;
ip.KeepUnmatched=true;
ip.addRequired('ND2filenameORreader')
ip.addParameter('numImagePairsToTryToFindAngleFor',20,@isnumeric)
ip.addParameter('tileIDpairs',[],@(x) and(isnumeric(x),size(x,2)==2))
ip.parse(ND2filenameORreader,varargin{:});

ND2filenameORreader=ip.Results.ND2filenameORreader;
numImagePairsToTryToFindAngleFor=ip.Results.numImagePairsToTryToFindAngleFor;
tileIDpairs=ip.Results.tileIDpairs;

% handle empty varargin case
if isempty(varargin)
    % skip input parsing
    vararginPassthrough={};
else %something was in varargin, pass on what wasn't used
    varginPassthrough=ip.Unmatched;
end



%% tileIDlist
numTilesInNd2file=size(stagecoor_XY,1);
tileIDlist=[1:numTilesInNd2file]';

%% numImagePairsToTryToFindAngleFor
if numImagePairsToTryToFindAngleFor>numTilesInNd2file
    numImagePairsToTryToFindAngleFor=numTilesInNd2file-1;
end

%% What tileIDpairs to look at
%tileIDlist=1:numImagePairsToTryToFindAngleFor;
if isempty(tileIDpairs)
    % generate a random tileIDpairs matrix
    rng default;
    idxChosen=randperm(length(tileIDlist)-1,numImagePairsToTryToFindAngleFor)';
    tile1IDsChosen=tileIDlist(idxChosen);
    tile2IDsChosen=tile1IDsChosen+1;
    tileIDpairs=[tile1IDsChosen, tile2IDsChosen];
elseif size(tileIDpairs,1)<numImagePairsToTryToFindAngleFor
    % pad the input with more random pairs
    tileIDsToChooseFrom=tileIDlist(~ismember(tileIDlist,tileIDpairs(:,1)));
    tileIDsToChooseFrom=tileIDsToChooseFrom(tileIDsToChooseFrom~=tileIDlist(end)); % not the last one since second image is +1;
    rng default;
    idxChosen=randperm(length(tileIDsToChooseFrom),numImagePairsToTryToFindAngleFor-size(tileIDpairs,1))';
    tile1IDsChosen=tileIDsToChooseFrom(idxChosen);
    tile2IDsChosen=tile1IDsChosen+1;
    tileIDpairs=[tileIDpairs;[tile1IDsChosen, tile2IDsChosen]];
elseif size(tileIDpairs,1)>numImagePairsToTryToFindAngleFor
    numImagePairsToTryToFindAngleFor=size(tileIDpairs,1);
else
    error('why am i here')
end

assert(numImagePairsToTryToFindAngleFor==size(tileIDpairs,1))

% loop through tileIDpairs
cameraAngleList=nan(numImagePairsToTryToFindAngleFor,1);
statusFlagList=nan(numImagePairsToTryToFindAngleFor,1);
numInliersList=nan(numImagePairsToTryToFindAngleFor,1);
fprintf('  findCameraAngle using %i image pairs:\n',numImagePairsToTryToFindAngleFor)
    
for iTileIDpair=1:numImagePairsToTryToFindAngleFor
    tileID1=tileIDpairs(iTileIDpair,1);
    tileID2=tileIDpairs(iTileIDpair,2);
    
    uint16_img1= getPlaneFromNd2file(ND2filenameORreader, tileID1, 'DAPI','closeReaderAfterGettingPlane',false);
    uint16_img2=getPlaneFromNd2file(ND2filenameORreader, tileID2, 'DAPI','closeReaderAfterGettingPlane',false);
    
    % convert uint16 arrays to matlab images
    img1=mat2gray(uint16_img1,[0 65535]);
    img2=mat2gray(uint16_img2,[0 65535]);
    
    % % add rotation/flipping to test this algorithm
    % img1_beforeTransform=img1;
    %img1=imrotate(img1,90);
    %img2=imrotate(img2,90);
    %img1=flip(img1,1); % 1=row dimension reversed, 2=col dimension reversed
    %img2=flip(img2,1); % 1=row dimension reversed, 2=col dimension reversed
    % imshowpair(img1_beforeTransform,img1,'Scaling','joint');
    % title('before (green), after (purple)')
    
    
    DisplayMatchingPointPairsUsedInTform=false; %shows figure of point pairs
    
    stagecoord_xy1_um=stagecoor_XY(tileID1,:);
    stagecoord_xy2_um=stagecoor_XY(tileID2,:);
    
    %[cameraAngle,imgcoord_tform,statusFlag]=findCameraAngleBetweenTwoImages(img1,img2,stagecoor_xy1_um,stage_xy2_um,UmPerPixel,StageAndImageAxesAreMirrored,'DisplayMatchingPointPairsUsedInTform',DisplayMatchingPointPairsUsedInTform,'MetricThreshold',20,'NumberOfStrongestSURFFeaturesToUse',1000);
    [cameraAngleThisTilePair,imgcoord_tform,numInliers,statusFlag]=findCameraAngleBetweenTwoSameSizeImages(img1,img2,stagecoord_xy1_um,stagecoord_xy2_um,UmPerPixel,StageAndImageAxesAreMirrored,varginPassthrough);
    % statusFlag from estimateGeometricTransform2D:
    %   0: No error.
    %   1: matchedPoints1 and matchedPoints2 do not contain enough points.
    %   2: Not enough inliers have been found
    % statusFlag from esimateTform:
    %   3: not enough inliers
    % statusFlag from findCameraAngleBetweenTwoImages:
    %   4: check whether StageYDimensionIsNotFlipped
    %
    % imgcoord_tform =
    %
    %       rigid2d with properties:
    %
    %            Rotation: [2×2 single]
    %         Translation: [1.8404e+03 -0.5713]
    

    fprintf('      from img pair %3i (tileIDs %4i and %4i): angle = %3.03f, %3i inliers, statusFlag=%i\n',iTileIDpair,tileID1,tileID2,cameraAngleThisTilePair,numInliers,statusFlag)

    cameraAngleList(iTileIDpair)=cameraAngleThisTilePair;
    numInliersList(iTileIDpair)=numInliers;
    statusFlagList(iTileIDpair)=statusFlag;
end % end tile loop

cameraAngleListNotnan=cameraAngleList(~isnan(cameraAngleList));
cameraAngle=median(cameraAngleListNotnan);

    fprintf('    final result: cameraAngle=%.03f\n',cameraAngle)
end