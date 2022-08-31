function T=getUprightFullTileImageTable(tileIDs,numImgRows,numImgCols,cameraAngleVsRefAngle)
% simplified version of getFullTileImageTable where image coordinates
% always reflect the default directionalities of rscan. Ie. image
% flipping and rotation must be handled right when grabbing the image tile,
% then after those transformations these coordinates apply.

numTiles=size(tileIDs,1);

cameraAngleVsRefAngle=wrapTo360(cameraAngleVsRefAngle);


% if rowDimIsFlippedVsYDim_refRound==rowDimIsFlippedVsYDim_thisRound
%    RowOne=1;
%    RowEnd=numImgRows;
%    ColOne=1;
%    ColEnd=numImgCols;
% else
%    RowOne=numImgRows;
%    RowEnd=1;
%    ColOne=1;
%    ColEnd=numImgCols;    
% end
% 
% if any(cameraAngleVsRefAngle==[0 360])
%     TLrow=RowOne;
%     TLcol=ColOne;
%     BRrow=RowEnd;
%     BRcol=ColEnd;
% elseif cameraAngleVsRefAngle==90
%     TLrow=RowEnd;
%     TLcol=ColOne;
%     BRrow=RowOne;
%     BRcol=ColEnd;
% elseif cameraAngleVsRefAngle==180
%     TLrow=RowEnd;
%     TLcol=ColEnd;
%     BRrow=RowOne;
%     BRcol=ColEnd;
% elseif cameraAngleVsRefAngle==270
%     TLrow=RowOne;
%     TLcol=ColEnd;
%     BRrow=RowEnd;
%     BRcol=ColOne;
% else
%     error('error in logic')
% end
% 
% 
% T=       array2table([tileIDs ,repmat(TLrow,numTiles,1),repmat(TLcol,numTiles,1),repmat(BRrow,numTiles,1),repmat(BRcol,numTiles,1)],...
%     'VariableNames',{'tileID',       'TLrow',                 'TLcol',                 'BRrow',                 'BRcol'});

if any(cameraAngleVsRefAngle==[90 270]) % switch dimensions
    temp=numImgRows;
    numImgRows=numImgCols;
    numImgCols=temp;
end

T=      array2table([tileIDs,     repmat(1,numTiles,1),     repmat(numImgRows,numTiles,1), repmat(1,numTiles,1),     repmat(numImgCols,numTiles,1)       ],...
    'VariableNames',{'tileID',         'TopRow',                 'BottomRow',                     'LeftCol',              'RightCol'});

end