function [cameraAngle,imgcoord_tform,numInliers,statusFlag]=findCameraAngleBetweenTwoSameSizeImages(img1,img2,stagecoord_xy1_um,stagecoord_xy2_um,UmPerPixel,rowDimIsFlippedVsYDim,varargin)
% [cameraAngle,imgcoord_tform,statusFlag]=findCameraAngleBetweenTwoImages(img1,img2,stagecoord_xy1_um,stagecoord_xy2_um,UmPerPixel,StageAndImageAxesAreMirrored,varargin)
% 
% findCameraAngleBetweenTwoImages
%
% can call this with more high-level findCameraOrientation (formerly findCameraAngle) with ND2 file input
% 
% vargins get passed through to estimateTform
% type the command 'help estimateTform' to see what optional arguments are
% possible
% 
% StageAndImageAxesAreMirrored: 
%       If you force the stagecoord +X dimension to point in the
%       same direction as the imagecoord +col dimension, then:
%            if  the stagecoord +Y dimension points in the +row dimension,
%                   then StageAndImageAxesAreMirrored=true (should be default)
%            if the stagecoord +Y dimension points in the -row dimension,
%                   then StageAndImageAxesAreMirrored=false
%
% statusFlag from estimateGeometricTransform2D:
%   0: No error.
%   1: matchedPoints1 and matchedPoints2 do not contain enough points.
%   2: Not enough inliers have been found (per estimateGeometricTransform2D)
% statusFlag from this function:
%   3: not enough inliers (per user-supplied or default minNumberOfInliers)
%

assert(isequal(size(img1),size(img2)))

ip=inputParser;
ip.KeepUnmatched=true;
ip.addParameter('minNumberOfInliers',2,@isnumeric)

ip.parse(varargin{:})
minNumberOfInliers=ip.Results.minNumberOfInliers;

if isempty(varargin)
    % skip input parsing
    varginPassthrough=struct();
else %something was in varargin, pass on what wasn't used
    varginPassthrough=ip.Unmatched; % struct
end

% default values of 


% ip=inputParser;
% ip.addParameter('MetricThreshold',30);
%
% % handle empty varargin case
% for i=1:numel(varargin)
%     if  isempty(varargin{i})
%         keepIdx=true(1,numel(varargin)); keepIdx(i)=false; varargin=varargin(keepIdx);
%     end
% end

% are images even theoretically close enough to eachother
stageDistanceBetweenImages=sqrt(sum((stagecoord_xy1_um - stagecoord_xy2_um).^2));
maxAllowableDistOfImagesInUm=UmPerPixel*(max(size(img1))/2 + max(size(img2)));

if stageDistanceBetweenImages>maxAllowableDistOfImagesInUm
    cameraAngle=nan;
    imgcoord_tform=nan;
    statusFlag=4;
    return
end


%% call estimateTform
[imgcoord_tform,numInliers,statusFlag]=estimateTform(img2,img1,varginPassthrough); %(moveImg, refImg,...). imgcoord_tform transforms from moveImg-->fixedImg

%% enough inliers for a reliable transform?
if (statusFlag==0)&&(numInliers<minNumberOfInliers)
    statusFlag=3; % not enough inliers
end


if nargout<3 % user is not asking for statusFlag, alert them
    if statusFlag==1
        warning('statusFlag=1: matchedPoints1 and matchedPoints2 do not contain enough points.')
    elseif statusFlag==2
        warning('statusFlag=2: Not enough inliers have been found (per estimateGeometricTransform2D');
    elseif statusFlag==3
        warning("statusFlag=2: Not enough inliers have been found (per findCameraAngleBetweenTwoImages's minimum of %i",minNumberOfInliers);
    end
end

if statusFlag>0
    cameraAngle=nan;
    imgcoord_tform=nan;
    return
end



%% convert the imgcoord_tform and stage coordinates into a camera angle
imgcoord_translation_Col=double(imgcoord_tform.Translation(1)); % add this to move X to get fixed X
imgcoord_translation_Row=double(imgcoord_tform.Translation(2)); % add this to move Y to get fixed Y
imgcoord_translation_Magnitude=norm([imgcoord_translation_Col,imgcoord_translation_Row]);

% stage coordinates
stagecoord_movement_x=stagecoord_xy2_um(1)-stagecoord_xy1_um(1);
stagecoord_movement_y=stagecoord_xy2_um(2)-stagecoord_xy1_um(2);


% check that the magnitude of the translation is about what we expect
% stage coordinates --> Convert to expected imgcoord translation assuming angle=zero
% actually since all we use is magnitude here don't really have to assume zero angle
expected_imgcoord_translation_Col_ifAngleIsZero=(stagecoord_movement_x)/UmPerPixel;
expected_imgcoord_translation_Row_ifAngleIsZero=(stagecoord_movement_y)/UmPerPixel;
expected_imgcoord_translation_Magnitude_ifAngleIsZero=norm([expected_imgcoord_translation_Col_ifAngleIsZero,expected_imgcoord_translation_Row_ifAngleIsZero]);
ToleranceOfMagnitudeOfTranslationError=0.05; % 0.05 is 5% tolerance.
if ~and(imgcoord_translation_Magnitude*(1-ToleranceOfMagnitudeOfTranslationError)<=expected_imgcoord_translation_Magnitude_ifAngleIsZero,...
        imgcoord_translation_Magnitude*(1+ToleranceOfMagnitudeOfTranslationError)>=expected_imgcoord_translation_Magnitude_ifAngleIsZero)
    % transform from features is too far off from stage movement to be real
    statusFlag=4;
    cameraAngle=nan;
    imgcoord_tform=nan;
    return
end


% angle of camera implied by imgcoord translation
% angle=atan2(Y,X)*180/pi()
if rowDimIsFlippedVsYDim
    angleDeg_from_imgcoord_translation=atan2( -1* imgcoord_translation_Row,     imgcoord_translation_Col)*180/pi();
else
    angleDeg_from_imgcoord_translation=atan2( -1* imgcoord_translation_Row, -1* imgcoord_translation_Col)*180/pi();
end

% angle of stage coordinates from img1 --> img2
angleDeg_of_stagecoord_movement=atan2(stagecoord_movement_y,stagecoord_movement_x)*180/pi();

% cameraAngle = difference of these angles
cameraAngle=angdiff(pi()/180*angleDeg_from_imgcoord_translation,pi()/180*angleDeg_of_stagecoord_movement)*180/pi();
% example:
% angleDifference=angdiff(80*pi()/180,90*pi()/180)*180/pi()
% angleDifference = 10
end