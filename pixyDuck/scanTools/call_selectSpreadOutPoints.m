%call_selectedSpreadOutPoints

XYall=Scan.Rounds(1).fullTiles.stage.T{:,{'centerX','centerY'}};
indSelected=[Scan.Rounds(1).ControlPoints.tileID_move]'
numAdditional=5
ind=selectSpreadOutPoints(XYall,numAdditional,indSelected)
size(ind)