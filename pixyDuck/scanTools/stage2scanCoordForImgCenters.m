function [regscancoord_RowOffsetFromREF,regscancoord_ColOffsetFromREF]= stage2scanCoordForImgCenters(regstagecoord_CenterX,regstagecoord_CenterY,REF_regstagecoord_ULtile_CenterX,REF_regstagecoord_ULtile_CenterY,cameraAngle,UmPerPixel,StageAndImageAxesAreMirrored)
% [regscancoord_RowOffsetFromREF,regscancoord_ColOffsetFromREF]= stage2scanCoordForImgCenters(regstagecoord_CenterX,regstagecoord_CenterY,REF_regstagecoord_ULtile_CenterX,REF_regstagecoord_ULtile_CenterY,cameraAngle,UmPerPixel)
%

% this function definitely needs work

% add two components to the offset:
%   component 1 = from natural stage XY offset
%   component 2 = from camera angle not being on-axis
% Not an amazing way to do it but...

assert(StageAndImageAxesAreMirrored==true) % not sure how it works otherwise
assert(and(cameraAngle<10,cameraAngle>-10)) % not sure how it works otherwise

% component 1 from stage offsets
regstagecoord_OffsetXfromREF=regstagecoord_CenterX - REF_regstagecoord_ULtile_CenterX;
regstagecoord_OffsetYfromREF=regstagecoord_CenterY - REF_regstagecoord_ULtile_CenterY;


component1_regscancoord_RowOffsetFromRef =   regstagecoord_OffsetYfromREF/UmPerPixel;
component1_regscancoord_ColOffsetFromRef = - regstagecoord_OffsetXfromREF/UmPerPixel;

% component 2 from camera angle
regstagecoord_HypotenusefromREF=sqrt(regstagecoord_OffsetXfromREF.^2 + regstagecoord_OffsetYfromREF.^2);

component2_regscancoord_RowOffsetFromREF =  sin(cameraAngle*pi()/180)      * regstagecoord_HypotenusefromREF/UmPerPixel;
component2_regscancoord_ColOffsetFromREF = (cos(cameraAngle*pi()/180) - 1) * regstagecoord_HypotenusefromREF/UmPerPixel;


regscancoord_RowOffsetFromREF= component1_regscancoord_RowOffsetFromRef + component2_regscancoord_RowOffsetFromREF;
regscancoord_ColOffsetFromREF= component1_regscancoord_ColOffsetFromRef + component2_regscancoord_ColOffsetFromREF;

end