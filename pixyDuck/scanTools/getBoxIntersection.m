function TboxIntersect=getBoxIntersection(T)
%         T=array2table([tileIDs   TopRow   BottomRow   LeftCol   RightCol],...
%      'VariableNames',{'tileIDs','TopRow','BottomRow','LeftCol','RightCol'});

% intersection
TopRow=max(T.TopRow);
BottomRow=min(T.BottomRow);

LeftCol=max(T.LeftCol);
RightCol=min(T.RightCol);

TboxIntersect=array2table([TopRow   BottomRow   LeftCol   RightCol],...
         'VariableNames',{'TopRow','BottomRow','LeftCol','RightCol'});

end