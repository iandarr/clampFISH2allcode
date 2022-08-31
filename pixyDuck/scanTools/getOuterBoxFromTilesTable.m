function TouterBox=getOuterBoxFromTilesTable(T)

%T=          array2table( [tileIDs, CenterRowExact,  CenterColExact,  CenterRow,  CenterCol,  TopRow,  BottomRow,  LeftCol,  RightCol,  numScanRows,  numScanCols],...
     %    'VariableNames',{'TileID','CenterRowExact','CenterColExact','CenterRow','CenterCol','TopRow','BottomRow','LeftCol','RightCol','numScanRows','numScanCols'});


TopRow=min(T.TopRow);
BottomRow=max(T.BottomRow);
LeftCol=min(T.LeftCol);
RightCol=max(T.RightCol);


TouterBox=array2table([ TopRow   BottomRow   LeftCol   RightCol],...
     'VariableNames', {'TopRow','BottomRow','LeftCol','RightCol'});
end