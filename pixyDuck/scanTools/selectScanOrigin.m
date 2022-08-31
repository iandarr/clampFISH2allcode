function originXY=selectScanOrigin(T,cameraAngle,rowDimIsFlippedVsYDim)
% originXY=selectScanOrigin(Scan.Rounds(referenceRound).fullTiles.rstage.T)

cameraAngle=wrapTo360(cameraAngle);

XYall=[T{:,{'X_RowOneColOne','Y_RowOneColOne'}};...
       T{:,{'X_RowEndColOne','Y_RowEndColOne'}};...
       T{:,{'X_RowOneColEnd','Y_RowOneColEnd'}};...
       T{:,{'X_RowEndColEnd','Y_RowEndColEnd'}}];

R = [cosd(cameraAngle) -sind(cameraAngle); sind(cameraAngle) cosd(cameraAngle)];

%% projection with fzero way
colDir=R*[1;0];
if rowDimIsFlippedVsYDim
    rowDir=R*[0;-1];
else
    rowDir=R*[0;1];
end

[~,idxColDirMin]=min(XYall*colDir);
[~,idxRowDirMin]=min(XYall*rowDir);
%originXY=[XYall(idxColDirMin,1),XYall(idxRowDirMin,2)];
xyColDirMin=XYall(idxColDirMin,:);
xyRowDirMin=XYall(idxRowDirMin,:);

if any(cameraAngle==[90 270])
    % don't project
    originXY=[xyRowDirMin(1) xyColDirMin(2)];
elseif any(cameraAngle==[0 180])
    originXY=[xyColDirMin(1) xyRowDirMin(2)];
else
    % project these and get intersection point

p1=polyfit([xyColDirMin(1),xyColDirMin(1)+rowDir(1)],[xyColDirMin(2),xyColDirMin(2)+rowDir(2)],1);
p2=polyfit([xyRowDirMin(1),xyRowDirMin(1)+colDir(1)],[xyRowDirMin(2),xyRowDirMin(2)+colDir(2)],1);
originX = fzero(@(x) polyval(p1-p2,x),xyColDirMin(1));
originY= polyval(p1,originX);
originXY=[originX,originY];

end
%% bounding rectangle way
% [rectxBound,rectyBound,~]=minboundrect(XYall(:,1),XYall(:,2)); % get minimal bounding rectangle
% 
% rectxy4Bound=[rectxBound(1:4),rectyBound(1:4)];
% 
% if rowDimIsFlippedVsYDim
%     biggerColRowDir=R*[1;-1];
% else
%     biggerColRowDir=R*[1; 1];
% end
% [~,idxRectxy4Bound]=min(rectxy4Bound*biggerColRowDir);
% originXY=rectxy4Bound(idxRectxy4Bound,:);



end