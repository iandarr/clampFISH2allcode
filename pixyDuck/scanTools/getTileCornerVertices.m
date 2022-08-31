function Tvertices=getTileCornerVertices(centerXY,numRows,numCols,UmPerPixel,cameraAngle,rowDimIsFlippedVsYDim)
% centerXY is a Nx2 array
%
% output
% Tvertices=[TverticesX,TverticesY];
% TverticesX is a Nx4 table, with the variables X_RowOneColOne, X_RowOneColEnd, X_RowEndColOne, X_RowEndColEnd
% TverticesY is a Nx4 table, with the variables Y_RowOneColOne, Y_RowOneColEnd, Y_RowEndColOne, Y_RowEndColEnd
% 

assert(and(numRows>=2, numCols>=2))

% get corner points of an image centered about 0,0;
% assume cameraAngle=0 for now, where +column dimension == +x dimension
% if rowDimIsFlippedVsYDim=true, RowOneColOne is at top left (small x, large y)
% if rowDimIsFlippedVsYDim=false,RowOneColOne is at bot left (small x, small y)

verticesXY4aroundOrigin=getVerticesXY4aroundOrigin(numRows,numCols,UmPerPixel,cameraAngle,rowDimIsFlippedVsYDim); %2x4

Xvertices=centerXY(:,1) + verticesXY4aroundOrigin(1,:); % Nx1 + 1x4 matlab knows what to do
Yvertices=centerXY(:,2) + verticesXY4aroundOrigin(2,:); % Nx1 + 1x4 matlab knows what to do

TverticesX=array2table(Xvertices,'VariableNames',{'X_RowOneColOne','X_RowOneColEnd','X_RowEndColOne','X_RowEndColEnd'});
TverticesY=array2table(Yvertices,'VariableNames',{'Y_RowOneColOne','Y_RowOneColEnd','Y_RowEndColOne','Y_RowEndColEnd'});
Tvertices=[TverticesX,TverticesY];
end

function verticesXY4aroundOrigin=getVerticesXY4aroundOrigin(numRows,numCols,UmPerPixel,cameraAngle,rowDimIsFlippedVsYDim)
% returns a 2x4 matrix
halfHeight=numRows/2 -0.5; % -0.5 becaues this is to center of the corner pixel
halfWidth=numCols/2 - 0.5; % -0.5 because this is to center of the corner pixel

% [x,y] format, in pixels
if rowDimIsFlippedVsYDim
    RowOneColOne=[ -halfWidth;  halfHeight]; %(small x, large y)
    RowOneColEnd=[  halfWidth;  halfHeight]; %(large x, large y)
    RowEndColOne=[ -halfWidth; -halfHeight]; %(small x, small y)
    RowEndColEnd=[  halfWidth; -halfHeight]; %(large x, small y)
else
    RowOneColOne=[ -halfWidth; -halfHeight]; %(small x, small y)
    RowOneColEnd=[  halfWidth; -halfHeight]; %(large x, small y)
    RowEndColOne=[ -halfWidth;  halfHeight]; %(small x, large y)
    RowEndColEnd=[  halfWidth;  halfHeight]; %(large x, large y)
end
verticesUnrotatedUnscaled=[RowOneColOne,RowOneColEnd,RowEndColOne,RowEndColEnd];

% scale by UmPerPixel
verticesUnscaled=UmPerPixel*verticesUnrotatedUnscaled;

% Rotate with rotation matrix (for rotating points around 0,0)
R = [cosd(cameraAngle) -sind(cameraAngle); sind(cameraAngle) cosd(cameraAngle)];

verticesXY4aroundOrigin=R*verticesUnscaled;
end