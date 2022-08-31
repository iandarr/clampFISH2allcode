% call selectScanOrigin


cameraAngle=Scan.cameraAngleOfReferenceRound;
cameraAngle=1
rscanOriginXY=selectScanOrigin(Scan.Rounds(referenceRound).fullTiles.rstage.T,cameraAngle,Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim)
Scan.rscanOriginXY=rscanOriginXY;
