function [cameraAngle,rowDimIsFlippedVsYDim,cameraAngleList,numInliersList,statusFlagList]=findCameraOrientation(ND2filenameORreader,XY,UmPerPixel,ImgRowsCols,varargin)
% cameraAngle=findCameraAngle(ND2filenameORreader,XY,UmPerPixel,StageYDimensionIsNotFlipped,estimateTformOptionalNameValueArguments)
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
%   rowDimIsFlippedVsYDim:
%       If you force the stagecoord +X dimension to point in the
%       same direction as the imagecoord +col dimension, then:
%            if  the stagecoord +Y dimension points in the +row dimension,
%                   then rowDimIsFlippedVsYDim=true (should be default)
%            if the stagecoord +Y dimension points in the -row dimension,
%                   then rowDimIsFlippedVsYDim=false
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
%               will try to look at 30 image pairs and
%               output the median angle. Ignored if tileIDpairs is used
%
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
ip.addParameter('numTilePairsToTryToFindAngleFor',30,@isnumeric)
ip.addParameter('numSuccessfulTilePairsAfterWhichToStopFindingAngle',8,@isnumeric)
ip.addParameter('tileIDpairsForAngleUserInput',[],@(x) and(isnumeric(x),size(x,2)==2))
ip.parse(ND2filenameORreader,varargin{:});

ND2filenameORreader=ip.Results.ND2filenameORreader;
numTilePairsToTryToFindAngleFor=ip.Results.numTilePairsToTryToFindAngleFor;
tileIDpairsForAngleUserInput=ip.Results.tileIDpairsForAngleUserInput;
numSuccessfulTilePairsAfterWhichToStopFindingAngle=ip.Results.numSuccessfulTilePairsAfterWhichToStopFindingAngle;

numTilePairsToTryToFindFlippedDimFor=20; % hard coded
numSuccessfulTilePairsAfterWhichToStopFindingAngleForCameraFlip=3; % hard coded

% handle empty varargin case
if isempty(varargin)
    % skip input parsing
    varginPassthrough=struct();
else %something was in varargin, pass on what wasn't used
    varginPassthrough=ip.Unmatched;
end

%% iTileList
numTilesTotal=size(XY,1);
iTileList=[1:numTilesTotal]';
%% What tile IDs and scan direction ('x' or 'y') to use for angle determination
% image pairs oriented in the other direction will then be used to determine rowDimIsFlippedVsYDim:
% if the angles found from that direction are different by ~180 degrees, then need to set rowDimIsFlippedVsYDim=true
[XDimConsecPairs,YDimConsecPairs,XDimNonConsecPairs,YDimNonConsecPairs]=getEstimatedOverlappingTilePairs(XY,UmPerPixel*min(ImgRowsCols)); % outputs 'x' or 'y'
%% What tileIDpairs to look at for angle determination and for row flip
consecutivePairBiasFactor=1.2;
[tileIDpairsForAngle,tileIDpairsForFlip,dimForAngleDetermination]=chooseTilePairsForAngleAndRowFlip(XDimConsecPairs,YDimConsecPairs,XDimNonConsecPairs,YDimNonConsecPairs,numTilePairsToTryToFindAngleFor,numTilePairsToTryToFindFlippedDimFor,consecutivePairBiasFactor);
% dimForAngleDetermination ('x' or 'y') is the 'primary' dimension. The tileIDpairsForAngle are oriented towards one another in this dimension

% warning('tileIDpairsForFlip manual override for testing')
% tileIDpairsForFlip=[15 16; 4 7; 7 14;tileIDpairsForFlip];

%% Add in user input tile pairs
tileIDpairsForAngle=addUserInputTilePairsToAutoChosenTilePairs(tileIDpairsForAngleUserInput,tileIDpairsForAngle,numTilePairsToTryToFindAngleFor);

% %% reversing these for testing purposes only
% warning('reversed the angle and flip pairs')
% temp=tileIDpairsForAngle;
% tileIDpairsForAngle=tileIDpairsForFlip;
% tileIDpairsForFlip=temp;
% dimForAngleDetermination='y';
numTilePairsToTryToFindAngleFor     =size(tileIDpairsForAngle,1); % may be smaller if don't have enough pairs
numTilePairsToTryToFindFlippedDimFor=size(tileIDpairsForFlip,1); % may be smaller if don't have enough pairs

%% simulation parameters
simulatedImagesOn=false;
simulatedExtraCameraAngle=89+179.771; %-179.771 real (thought before it was 179.7508 real)
simulatedRowFlip=false;
if simulatedImagesOn
    warning('SIMULATION is on: simulatedExtraCameraAngle=%.03f, simulatedRowFlip=%i',simulatedExtraCameraAngle,simulatedRowFlip)
end
%% find camera angle: loop through tileIDpairsForAngle
cameraAngleList=nan(numTilePairsToTryToFindAngleFor,1);
statusFlagList=nan(numTilePairsToTryToFindAngleFor,1);
numInliersList=nan(numTilePairsToTryToFindAngleFor,1);
fprintf('  findCameraOrientation: finding camera angle. Will stop after %i successful pairs or after having tried %i image pairs:\n',numSuccessfulTilePairsAfterWhichToStopFindingAngle,numTilePairsToTryToFindAngleFor)

rowDimIsFlippedVsYDimAssumption=true; % in my images, rowDimIsFlippedVsYDim=true
fprintf('   using rowDimIsFlippedVsYDimAssumption =%i\n',rowDimIsFlippedVsYDimAssumption)

for iTileIDpair=1:numTilePairsToTryToFindAngleFor
    tileID1=tileIDpairsForAngle(iTileIDpair,1);
    tileID2=tileIDpairsForAngle(iTileIDpair,2);
    
    uint16_img1= getPlaneFromNd2file(ND2filenameORreader, tileID1, 'DAPI','closeReaderAfterGettingPlane',false);
    uint16_img2=getPlaneFromNd2file(ND2filenameORreader, tileID2, 'DAPI','closeReaderAfterGettingPlane',false);
    
    % convert uint16 arrays to matlab images
    img1=mat2gray(uint16_img1,[0 65535]);
    img2=mat2gray(uint16_img2,[0 65535]);
    
    % SIMULATION of rotation and simulatedRowDimIsFlippedVsYDim
    % % add rotation/flipping to test this algorithm
    if simulatedImagesOn
        %img1_beforeTransform=img1;
        if (simulatedRowFlip==true) % my images are flipped, so in reality rowDimIsFlippedVsYDim=true
            %img1_beforeTransform=flip(img1_beforeTransform,1);
            img1=flip(img1,1); % 1=row dimension reversed, 2=col dimension reversed
            img2=flip(img2,1); % 1=row dimension reversed, 2=col dimension reversed
        end
        img1=imrotate(img1,-simulatedExtraCameraAngle,'crop');
        img2=imrotate(img2,-simulatedExtraCameraAngle,'crop');
%         imshowpair(img1_beforeTransform,img1,'Scaling','joint');
%         title('before (green), after (purple)')
    end
    
    xy1_um=XY(tileID1,:);
    xy2_um=XY(tileID2,:);
    
    if ~isfield(varginPassthrough,'intensityNormType')
        varginPassthrough.intensityNormType='SameScaling';
    end
    if ~isfield(varginPassthrough,'MetricThreshold')
        varginPassthrough.MetricThreshold=100;
    end
    
    %[cameraAngle,imgcoord_tform,statusFlag]=findCameraAngleBetweenTwoImages(img1,img2,stagecoor_xy1_um,stage_xy2_um,UmPerPixel,StageAndImageAxesAreMirrored,'DisplayMatchingPointPairsUsedInTform',DisplayMatchingPointPairsUsedInTform,'MetricThreshold',20,'NumberOfStrongestSURFFeaturesToUse',1000);
    [cameraAngleThisTilePair,imgcoord_tform,numInliers,statusFlag]=findCameraAngleBetweenTwoSameSizeImages(img1,img2,xy1_um,xy2_um,UmPerPixel,rowDimIsFlippedVsYDimAssumption,varginPassthrough);

    
    
    fprintf('      from img pair %3i (tileIDs %4i and %4i): angle = %3.03f, %3i inliers, statusFlag=%i\n',iTileIDpair,tileID1,tileID2,cameraAngleThisTilePair,numInliers,statusFlag)
    
    cameraAngleList(iTileIDpair)=cameraAngleThisTilePair;
    numInliersList(iTileIDpair)=numInliers;
    statusFlagList(iTileIDpair)=statusFlag;
    
    if sum(~(isnan(cameraAngleList)))==numSuccessfulTilePairsAfterWhichToStopFindingAngle
        fprintf('    since successfully found angle for %i pairs, stop finding angles for pairs\n',numSuccessfulTilePairsAfterWhichToStopFindingAngle);
        break
    end
end % end tile loop

%cameraAngleListNotnan=cameraAngleList(~isnan(cameraAngleList));
cameraAngle=middleAngle(cameraAngleList);
cameraAngle=wrapTo180(cameraAngle);

fprintf('    final result: cameraAngle=%.03f\n',cameraAngle)

%% determine StageAndImageAxesAreMirroredAssumption
% loop through tileIDpairsForFlip and also find angle. If this cameraAngle is is ~90 or ~270 degrees from the above cameraAngle, then assumed StageAndImageAxesAreMirroredAssumption is opposite
fprintf('  findCameraOrientation: checking StageAndImageAxesAreMirrored assumption. Will stop after %i successful pairs or after having tried %i image pairs:\n',numSuccessfulTilePairsAfterWhichToStopFindingAngleForCameraFlip,numTilePairsToTryToFindAngleFor)
for iTileIDpair=1:numTilePairsToTryToFindFlippedDimFor
    tileID1=tileIDpairsForFlip(iTileIDpair,1);
    tileID2=tileIDpairsForFlip(iTileIDpair,2);
    
    uint16_img1= getPlaneFromNd2file(ND2filenameORreader, tileID1, 'DAPI','closeReaderAfterGettingPlane',false);
    uint16_img2=getPlaneFromNd2file(ND2filenameORreader, tileID2, 'DAPI','closeReaderAfterGettingPlane',false);
    
    % convert uint16 arrays to matlab images
    img1=mat2gray(uint16_img1,[0 65535]);
    img2=mat2gray(uint16_img2,[0 65535]);
    
    % SIMULATION of rotation and simulatedRowDimIsFlippedVsYDim
    % % add rotation/flipping to test this algorithm
    if simulatedImagesOn
        %img1_beforeTransform=img1;
        if (simulatedRowFlip==true) % my images are flipped, so in reality rowDimIsFlippedVsYDim=true
            %img1_beforeTransform=flip(img1_beforeTransform,1);
            img1=flip(img1,1); % 1=row dimension reversed, 2=col dimension reversed
            img2=flip(img2,1); % 1=row dimension reversed, 2=col dimension reversed
        end
        img1=imrotate(img1,-simulatedExtraCameraAngle,'crop');
        img2=imrotate(img2,-simulatedExtraCameraAngle,'crop');
        %imshowpair(img1_beforeTransform,img1,'Scaling','joint');
        %title('before (green), after (purple)')
    end
    
    xy1_um=XY(tileID1,:);
    xy2_um=XY(tileID2,:);
    [cameraAngleThisTilePair,imgcoord_tform,numInliers,statusFlag]=findCameraAngleBetweenTwoSameSizeImages(img1,img2,xy1_um,xy2_um,UmPerPixel,rowDimIsFlippedVsYDimAssumption,varginPassthrough);
    fprintf('      from img pair %3i (tileIDs %4i and %4i): angle = %3.03f, %3i inliers, statusFlag=%i\n',iTileIDpair,tileID1,tileID2,cameraAngleThisTilePair,numInliers,statusFlag)
    
    cameraAngleListForFlip(iTileIDpair)=cameraAngleThisTilePair;
    numInliersListForFlip(iTileIDpair)=numInliers;
    statusFlagListForFlip(iTileIDpair)=statusFlag;
    
    if sum(~(isnan(cameraAngleListForFlip)))==numSuccessfulTilePairsAfterWhichToStopFindingAngleForCameraFlip
        fprintf('    since successfully found angle for %i pairs, stop finding angles for pairs\n',numSuccessfulTilePairsAfterWhichToStopFindingAngleForCameraFlip);
        break
    end
end

cameraAngleForFlip=middleAngle(cameraAngleListForFlip);
fprintf('    final result: cameraAngleForFlip=%.03f\n',cameraAngleForFlip)

%% determine if rowDimIsFlippedVsYDimAssumption is wrong and adjust cameraAngle and rowDimIsFlippedVsYDim accordingly
% when this is wrong, the cameraAngle found from primaryDim and
% secondaryDim tile pairings are ~180 degrees apart
angDiff=wrapTo180(180/pi()*angdiff(pi()/180*cameraAngleForFlip,pi()/180*cameraAngle));
absAngDiffFrom180=abs(180-abs(angDiff));
angTol=10; % hard coded
if absAngDiffFrom180<angTol
    rowDimIsFlippedVsYDim=~rowDimIsFlippedVsYDimAssumption; % it is opposite of assumption
    fprintf('  angular difference between x and y registrations is angDiff=%.04f, or %0.04f degrees away from 180 degrees. Therefore, assumption about rowDimIsFlippedVsYDim is wrong. Changing to rowDimIsFlippedVsYDim=%i\n',angDiff,absAngDiffFrom180,rowDimIsFlippedVsYDim)
    if strcmp(dimForAngleDetermination,'y')
        cameraAngle= wrapTo180(180+cameraAngle); % we also need to rotate the cameraAngle. But only if primaryDim was 'y'
        fprintf('  and since cameraAngle was based on tile pairs in the y dimension, need to rotate cameraAngle 180 degrees to %.03f\n',cameraAngle)
    else
        fprintf('  since cameraAngle was based on tile pairs in the x dimension, keeping cameraAngle of %.03f\n',cameraAngle)
    end
else
    rowDimIsFlippedVsYDim= rowDimIsFlippedVsYDimAssumption;
    fprintf('  angular difference between x and y registrations is angDiff=%.04f, which is close to 0 degrees. Therefore, assumption about rowDimIsFlippedVsYDim was correct. Keeping rowDimIsFlippedVsYDim=%i and keeping cameraAngle=%.03f\n',angDiff,rowDimIsFlippedVsYDim,cameraAngle)
end


end