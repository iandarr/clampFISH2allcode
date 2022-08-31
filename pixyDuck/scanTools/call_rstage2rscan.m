% call rstage2rscan
%         stageCenterXY=Scan.Rounds(thisRound).fullTiles.rstage.T{:,{'centerX','centerY'}};
 umPerPixel=Scan.Rounds(referenceRound).umPerPixel; % change to iRes
%         %           rstage2rscan(stageXY,rscanOriginXY,cameraAngle,rowDimIsFlippedVsYDim,UmPerPixel)
%         rscanRowCol=rstage2rscan(stageCenterXY,Scan.rscanOriginXY,Scan.cameraAngleOfReferenceRound,Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim,umPerPixel);
%         
iRes=1;
                stageCenterXY=Scan.Rounds(thisRound).fullTiles.rstage.T{:,{'centerX','centerY'}};
        %umPerPixel=Scan.resolutions(iRes);
        
        if rem(Scan.Rounds(thisRound).fullTileNumRows,2)==0
            numRowIsEven=true;
        else
            numRowIsEven=false;
        end
        if rem(Scan.Rounds(thisRound).fullTileNumCols,2)==0
            numColIsEven=true;
        else
            numColIsEven=false;
        end
        %           rstage2rscan(stageXY,rscanOriginXY,cameraAngle,rowDimIsFlippedVsYDim,UmPerPixel,numRowIsEven,numColIsEven)
        rscanRowCol=rstage2rscanCenters(stageCenterXY,Scan.rscanOriginXY,Scan.cameraAngleOfReferenceRound,Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim,umPerPixel,numRowIsEven,numColIsEven);
