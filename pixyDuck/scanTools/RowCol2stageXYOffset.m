function stageXYOffset=RowCol2stageXYOffset(RowColOffset,cameraAngle,rowDimIsFlippedVsYDim,UmPerPixel)
% MoveRowCol + RowColOffset = FixedRowCol
%
% let theta = negative camera Angle (radians)
% let alpha = angle of RowColOffset vector (radians)

MagPixels=norm(RowColOffset);


theta= -1* cameraAngle * pi/180;

if rowDimIsFlippedVsYDim
    alpha=atan2( -1* RowColOffset(1),RowColOffset(2));
else
    alpha=atan2(     RowColOffset(1),RowColOffset(2));
end



%alphadeg=alpha*180/pi
%thetadeg=theta*180/pi
%                                          deltaX              deltaY
stageXYOffset = UmPerPixel * MagPixels * [cos(theta-alpha), -1 * sin(theta-alpha)];


end