function Timage_overlap=cropImageCoords(Trscan_orig,Trscan_overlap,Timage_orig,nativePxPerRscanPx)


numImages=height(Trscan_overlap);

tileIDs=Trscan_overlap.tileID;

% all tileIDs in order?
assert(all(diff(Trscan_orig.tileID))>0) 
assert(all(diff(Trscan_overlap.tileID))>0) 
assert(all(diff(Timage_orig.tileID))>0)

% only the overlapping tiles
Trscan_orig=Trscan_orig(ismember(Trscan_orig.tileID,tileIDs),:);
Timage_orig=Timage_orig(ismember(Timage_orig.tileID,tileIDs),:);

% get offsets
ToffsetsRscanRes=table();
ToffsetsRscanRes.tileID=tileIDs;
ToffsetsRscanRes{:,{'TopRow','BottomRow','LeftCol','RightCol'}}=Trscan_overlap{:,{'TopRow','BottomRow','LeftCol','RightCol'}} - ...
                                                                   Trscan_orig{:,{'TopRow','BottomRow','LeftCol','RightCol'}};
                                                             
Timage_overlapExact=table();
Timage_overlapExact.tileID=Timage_orig.tileID;
Timage_overlapExact{:,{'TopRow','BottomRow','LeftCol','RightCol'}}=Timage_orig{:,{'TopRow','BottomRow','LeftCol','RightCol'}}+...
                                                              ToffsetsRscanRes{:,{'TopRow','BottomRow','LeftCol','RightCol'}}*nativePxPerRscanPx;

% round to nearest pixel. prioritize preserving scale over stitch at bottom/right edges

% first TopRow and LeftCol
Timage_overlap=table();
Timage_overlap.tileID=Timage_overlapExact.tileID;
Timage_overlap{:,{'TopRow','LeftCol'}}=round(Timage_overlapExact{:,{'TopRow','LeftCol'}});

Timage_roundError=table();
Timage_roundError.tileID=Timage_overlap.tileID;
Timage_roundError{:,{'TopRow','LeftCol'}}=Timage_overlap{:,{'TopRow','LeftCol'}} - Timage_overlapExact{:,{'TopRow','LeftCol'}};

% now add rounderror before getting BottomRow and RightCol
Timage_overlap{:,{'BottomRow','RightCol'}} = round(...
                                             Timage_overlapExact{:,{'BottomRow','RightCol'}} + ...
                                             Timage_roundError{:,     {'TopRow', 'LeftCol'}});
Timage_overlap=movevars(Timage_overlap,'BottomRow','After','TopRow');

Timage_overlap.numRows=Timage_overlap.BottomRow - Timage_overlap.TopRow + 1;
Timage_overlap.numCols=Timage_overlap.RightCol - Timage_overlap.LeftCol + 1;
end