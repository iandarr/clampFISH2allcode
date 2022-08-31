function rscanRowCol=rstage2rscan(stageXY,rscanOriginXY,cameraAngle,rowDimIsFlippedVsYDim,UmPerPixel)%,numRowIsEven,numColIsEven,isNativeRes)
% no rounding

% rotation matrix
R = [cosd(cameraAngle) -sind(cameraAngle); sind(cameraAngle) cosd(cameraAngle)];

colDir=R*[1;0];
if rowDimIsFlippedVsYDim
    rowDir=R*[0;-1];
else
    rowDir=R*[0;1];
end

XY=stageXY-rscanOriginXY;


rscanRowCol=[XY*rowDir  XY*colDir]/UmPerPixel;

% if isNativeRes
%     % then round the center point location at this stage
%     subtractRowCol=zeros(1,2);
%     if numRowIsEven
%         subtractRowCol(1)=0.5;
%     end
%     if numColIsEven
%         subtractRowCol(2)=0.5;
%     end
%     
%     rscanRowCol=round(rscanRowCol - subtractRowCol) + subtractRowCol;
%     
% end

end