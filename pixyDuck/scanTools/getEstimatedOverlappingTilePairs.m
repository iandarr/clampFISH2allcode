function [XDimConsecPairs,YDimConsecPairs,XDimNonConsecPairs,YDimNonConsecPairs]=getEstimatedOverlappingTilePairs(XY,smallestImageDimensionUm)
% estimates if tiles overlap based on their center XY coordinates and the
% smallest image dimension in Um
%
% since the camera angle is unknown, and since here we only have one image
% dimension, this is not truly a measure of whether they overlap. It would
% work exactly if the images were circles

numTiles=size(XY,1);
overlapsMat=pdist2(XY,XY,'euclidean')<smallestImageDimensionUm; % even if doesn't satisfy this it could still be overlapping at the corners. Or if one image dimension is much bigger. But this is good enough for these purposes
%overlapsMat=(1-diag(ones(1,numTiles))).*overlapsMat; % remove self overlaps
overlapsMat=(1-tril(ones(numTiles))).*overlapsMat; % remove self overlaps and duplicates

[tile1,tile2] = find(overlapsMat);
tileOverlapPairs=[tile1,tile2];
tileOverlapPairs=sort(tileOverlapPairs,2); % smaller tile on left
tileOverlapPairs=sortrows(tileOverlapPairs,1); % smaller tileIDs on top (left column)


% determine whether primarily a X or a Y direction overlap
XYdirOfPair=XY(tileOverlapPairs(:,2),:) - XY(tileOverlapPairs(:,1),:);
isXdimPair=abs(XYdirOfPair(:,1))>=abs(XYdirOfPair(:,2));

% determine if it's a consecutive pair (Ie. tile ID 3 and 4 touching means
% its consecutive
isConsecutivePair=1==abs(tileOverlapPairs(:,2)-tileOverlapPairs(:,1));

% would rather go with the consecutive pairs, since i trust the relative XY stage
% position to be more accurate with these. However, if 


XDimConsecPairs=tileOverlapPairs(all(   [ isXdimPair, isConsecutivePair],2),:);
YDimConsecPairs=tileOverlapPairs(all(   [~isXdimPair, isConsecutivePair],2),:);
XDimNonConsecPairs=tileOverlapPairs(all([ isXdimPair,~isConsecutivePair],2),:);
YDimNonConsecPairs=tileOverlapPairs(all([~isXdimPair,~isConsecutivePair],2),:);
end