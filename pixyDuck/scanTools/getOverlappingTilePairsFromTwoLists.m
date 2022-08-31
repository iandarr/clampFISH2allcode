function [Idx1,Idx2,areaOverlap]=getOverlappingTilePairsFromTwoLists(centerXY1,centerXY2,WidthHeight1,WidthHeight2,varargin)
% [Idx1,Idx2,areaOverlap]=getOverlappingTilePairsFromTwoLists(centerXY1,centerXY2,WidthHeight1,WidthHeight2)
%
% % [Idx1,Idx2,areaOverlap]=getOverlappingTilePairsFromTwoLists(centerXY1,centerXY2,WidthHeight1,WidthHeight2,minFractionalAreaOfSmallerTile)
%
% Eg: minFractionalAreaOfSmallerTile=0.6 means over 60% of the full tile
% area of the smaller of the tile lists must be overlapping to be output
%
% 
% height relates to Y, Width relates to X
%
%
% WidthHeight1=[1024 1022]

if ~isempty(varargin)
   minFractionalAreaOfSmallerTile=varargin{1};
   assert(isnumeric(minFractionalAreaOfSmallerTile));
   assert(minFractionalAreaOfSmallerTile>0)
else
    minFractionalAreaOfSmallerTile=-inf;
end

numTiles1=size(centerXY1,1);
numTiles2=size(centerXY2,1);


assert(isequal(size(WidthHeight1),[1 2]))
assert(isequal(size(WidthHeight2),[1 2]))

assert(all(WidthHeight1>0))
assert(all(WidthHeight2>0))

%overlapAbsDixtInX=abs(pdist2(centerXY1(:,1),centerXY2(:,1)));
%overlapAbsDixtInY=abs(pdist2(centerXY1(:,2),centerXY2(:,2)));

% rectangle is [x,y,width,height]
rects1=[centerXY1(:,1)-WidthHeight1(1)/2, centerXY1(:,2)-WidthHeight1(2)/2, repmat(WidthHeight1(1),numTiles1,1), repmat(WidthHeight1(2),numTiles1,1)];
rects2=[centerXY2(:,1)-WidthHeight2(1)/2, centerXY2(:,2)-WidthHeight2(2)/2, repmat(WidthHeight2(1),numTiles2,1), repmat(WidthHeight2(2),numTiles2,1)];
areasMat = rectint(rects1,rects2);
% areasMat is size numTiles1 x numTiles2

IdxMat=areasMat>0;
[Idx1,Idx2]=find(IdxMat);
areaOverlap=areasMat(IdxMat);

% take out indices where area of overlap is too small per minFractionalAreaOfSmallerTile
minArea=minFractionalAreaOfSmallerTile*min(WidthHeight1(1)*WidthHeight1(2),WidthHeight2(1)*WidthHeight2(2));
areaOfOverlapIsLargeEnough=areaOverlap>=minArea;
Idx1=Idx1(areaOfOverlapIsLargeEnough);
Idx2=Idx2(areaOfOverlapIsLargeEnough);
areaOverlap=areaOverlap(areaOfOverlapIsLargeEnough);
end