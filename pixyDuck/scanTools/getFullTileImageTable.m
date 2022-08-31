function T=getFullTileImageTable(tileIDs,numImgRows,numImgCols,cameraAngleVsRefAngle,rowDimIsFlippedVsYDim_refRound,rowDimIsFlippedVsYDim_thisRound)

numTiles=size(tileIDs,1);

cameraAngleVsRefAngle=wrapTo360(cameraAngleVsRefAngle);

if rowDimIsFlippedVsYDim_refRound==rowDimIsFlippedVsYDim_thisRound
   RowOne=1;
   RowEnd=numImgRows;
   ColOne=1;
   ColEnd=numImgCols;
else
   RowOne=numImgRows;
   RowEnd=1;
   ColOne=1;
   ColEnd=numImgCols;    
end

if any(cameraAngleVsRefAngle==[0 360])
    TLrow=RowOne;
    TLcol=ColOne;
    BRrow=RowEnd;
    BRcol=ColEnd;
elseif cameraAngleVsRefAngle==90
    TLrow=RowEnd;
    TLcol=ColOne;
    BRrow=RowOne;
    BRcol=ColEnd;
elseif cameraAngleVsRefAngle==180
    TLrow=RowEnd;
    TLcol=ColEnd;
    BRrow=RowOne;
    BRcol=ColEnd;
elseif cameraAngleVsRefAngle==270
    TLrow=RowOne;
    TLcol=ColEnd;
    BRrow=RowEnd;
    BRcol=ColOne;
else
    error('error in logic')
end

T=       array2table([tileIDs ,repmat(TLrow,numTiles,1),repmat(TLcol,numTiles,1),repmat(BRrow,numTiles,1),repmat(BRcol,numTiles,1)],...
    'VariableNames',{'tileID',       'TLrow',                 'TLcol',                 'BRrow',                 'BRcol'});

end