function scanMatrix=XYtoScanMatrix(XYpos,RowIsDim,ColIsDim)
% Scan.Rounds(round).stagecoordXYpos,Scan.RowIsDim,Scan.ColIsDim)
%
% XYpos is 2 column vector of tile stage positions (first column is X, second column is Y)
% RowIsDim='-y'
% ColIsDim='+x'

numTile=size(XYpos,1);

if numTile==1
    scanMatrix=1;
elseif numTile==0
    error('first input XYpos should be 2 columns with 1 tile each')
else % normal case with many tiles
    % look at first 2 tiles to get the expected regular inter-tile spacing (in stage coordinates)
    XYdist=XYpos(2,:)-XYpos(1,:);
    [gridSpacingGuess,scanDirectionDimension]=max(abs(XYdist));
    if scanDirectionDimension==1 % X dimension
        if strcmpi(RowIsDim,'x')
        scanDirectionGuessStageCoords='left'
        
    end
    
    for iTile=1:numTile
        
        
    end
    
end


end