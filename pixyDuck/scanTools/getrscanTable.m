function T=getrscanTable(tileIDs,rscanRowColCentersExact,numImgRows,numImgCols,nativePxPerRscanPx,cameraAngleVsRefAngle)

numTiles=size(tileIDs,1);

CenterRowExact=rscanRowColCentersExact(:,1);
CenterColExact=rscanRowColCentersExact(:,2);

assert(rem(cameraAngleVsRefAngle,90)==0);
assert(and(numImgRows>=1,numImgCols>=1))

% convert image to scan rowcols
if any([90 270]==wrapTo360(cameraAngleVsRefAngle))
    numScanRows=max(round(numImgCols/nativePxPerRscanPx),1);
    numScanCols=max(round(numImgRows/nativePxPerRscanPx),1);
else
    numScanRows=max(round(numImgRows/nativePxPerRscanPx),1);
    numScanCols=max(round(numImgCols/nativePxPerRscanPx),1);
end


% choose CenterRow & CenterCol (rounded to nearest pixel, but ends in 0.5 if even number of rows/cols)
if rem(numScanRows,2)==0
    % number of rows is even
    CenterRow=round(CenterRowExact - 0.5) + 0.5;
else
    CenterRow=round(CenterRowExact);
end

if rem(numScanCols,2)==0
    % number of cols is even
    CenterCol=round(CenterColExact - 0.5) + 0.5;
else
    CenterCol=round(CenterColExact);
end


TopRow=CenterRow - (numScanRows-1)/2;
BottomRow=TopRow + numScanRows -1;

LeftCol=CenterCol - (numScanCols-1)/2;
RightCol=LeftCol + numScanCols - 1;

T=      array2table( [tileIDs, CenterRowExact,  CenterColExact,  CenterRow,  CenterCol,  TopRow,  BottomRow,  LeftCol,  RightCol,  repmat(numScanRows,numTiles,1),  repmat(numScanCols,numTiles,1)],...
    'VariableNames',{'tileID','CenterRowExact','CenterColExact','CenterRow','CenterCol','TopRow','BottomRow','LeftCol','RightCol',       'numScanRows',                   'numScanCols'});

end