function [scancoord_offsetRow,scancoord_offsetCol]=cameraAngle2UpperLeftCornerOffets(cameraAngle,StageAndImageAxesAreMirrored,imgHeight,imgWidth)
% [scancoord_offsetRow,scancoord_offsetCol]=angleToCornerShift(cameraAngle,StageAndImageAxesAreMirrored,imgHeight,imgWidtht,imgCorner)
%
% cameraAngle in degrees
% must be upper left corner, not built for others right now 

assert(StageAndImageAxesAreMirrored)
% haven't worked out how it works otherwise


halfImgToMidpixelHeight=imgHeight/2 -0.5;
halfImgToMidpixelWidth=imgWidth/2 - 0.5;
hypotenuseToMidpixelCorner=sqrt(halfImgToMidpixelHeight^2 + halfImgToMidpixelWidth^2);

%% assuming no camera tilt, Row & Col:
% angle=atan2(Y,X)*180/pi()
angleToCorner_NoCameraTilt=atan2(halfImgToMidpixelHeight,-1*halfImgToMidpixelWidth)*180/pi();
% if we define center of image as scancoord row=0, col=0;
scancoord_UpperLeftRow_NoCameraTilt=-1*hypotenuseToMidpixelCorner*sin(angleToCorner_NoCameraTilt*pi()/180);
scancoord_UpperLeftCol_NoCameraTilt=   hypotenuseToMidpixelCorner*cos(angleToCorner_NoCameraTilt*pi()/180);

%% with camera tilt, what is angle
angleToCorner_WithCameraTilt=wrapTo180(angleToCorner_NoCameraTilt + cameraAngle);

%% with camera tilt, what are new stagecoords
scancoord_UpperLeftRow_WithCameraTilt=-1*hypotenuseToMidpixelCorner*sin(angleToCorner_WithCameraTilt*pi()/180);
scancoord_UpperLeftCol_WithCameraTilt=   hypotenuseToMidpixelCorner*cos(angleToCorner_WithCameraTilt*pi()/180);

%% what are the offsets (WithTilt - NoTilt)
scancoord_offsetRow= scancoord_UpperLeftRow_WithCameraTilt - scancoord_UpperLeftRow_NoCameraTilt;
scancoord_offsetCol= scancoord_UpperLeftCol_WithCameraTilt - scancoord_UpperLeftCol_NoCameraTilt;


end