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

assert(isequal(size(img1),size(img2)))

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



[imgcoord_tform,numInliers,statusFlag]=estimateTform(img1,img2,varargin{:}); %(moveImg, refImg,...)
if nargout<3 % user is not asking for statusFlag, alert them
if statusFlag==1
    warning('statusFlag=1: matchedPoints1 and matchedPoints2 do not contain enough points.')
elseif statusFlag==2
    warning('statusFlag=2: Not enough inliers have been found');
end
end

if statusFlag>0
    cameraAngle=nan;
    imgcoord_tform=nan;
    return
end
% statusFlag from estimateGeometricTransform2D:
%   0: No error.
%   1: matchedPoints1 and matchedPoints2 do not contain enough points.
%   2: Not enough inliers have been found
% statusFlag from esimateTform:
%   3: not enough inliers
%
% 
% imgcoord_tform =
%   rigid2d with properties:
%
%        Rotation: [2×2 single]
%     Translation: [1.8404e+03 -0.5713]
%  for example


%% convert the imgcoord_tform and stage coordinates into a camera angle
imgcoord_translation_Col=double(imgcoord_tform.Translation(1));
imgcoord_translation_Row=double(imgcoord_tform.Translation(2));
%imgcoord_translation_Euclidean=sqrt(imgcoord_translation_Col^2 + imgcoord_translation_Row^2);

% % which imgcoord dimention changed the most (col or row)
% if abs(imgcoord_translation_Col)>imgcoord_translation_Row
%     major_imcoord_translation_dim='col';
% else
%     major_imcoord_translation_dim='row';
% end


stagecoord_movement_x=stagecoord_xy2_um(1)-stagecoord_xy1_um(1);
stagecoord_movement_y=stagecoord_xy2_um(2)-stagecoord_xy1_um(2);

if rowDimIsFlippedVsYDim
    % then stage Y dimension is flipped, so need to negate stagecoord_movement_y
    stagecoord_movement_y= -1 * stagecoord_movement_y;   
end

% which stagecoord dimention changed the most (x or y)
% if abs(stagecoord_movement_x)>abs(stagecoord_movement_y)
%     major_stagecoord_movement_dim='x';
% else
%     major_stagecoord_movement_dim='y';
% end


% if strcmp(major_imcoord_translation_dim,'row') % images shifted primarily in row
%     if major_stagecoord_movement_dim==1 % scan direction is in x or -x
%         imgcoordRowIsApproximatelyStageDirection='+x';
%         imgcoordColIsApproximatelyStageDirection
%     else %major_stagecoord_movement_dim==2 % scan direction is in y or -y
%         
%     end
%     
% elseif strcmp(major_imcoord_translation_dim,'col') % images shifted primarily in row
%     
% else error('error in my logic')
% end


% angle of camera implied by imgcoord translation
% angle=atan2(Y,X)*180/pi()
angleDeg_from_imgcoord_translation=atan2(imgcoord_translation_Row,imgcoord_translation_Col)*180/pi();

% angle of stage coordinates from img1 --> img2
angleDeg_of_stagecoord_movement=atan2(stagecoord_movement_y,stagecoord_movement_x)*180/pi();
% example:
% angleDifference=angdiff(80*pi()/180,90*pi()/180)*180/pi()
% angleDifference = 10

cameraAngle=angdiff(pi()/180*angleDeg_of_stagecoord_movement,pi()/180*angleDeg_from_imgcoord_translation)*180/pi();

end