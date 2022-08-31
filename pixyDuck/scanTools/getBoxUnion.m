function TboxUnion=getBoxUnion(T)
%         T=array2table([tileIDs   TopRow   BottomRow   LeftCol   RightCol],...
%      'VariableNames',{'tileIDs','TopRow','BottomRow','LeftCol','RightCol'});

% union
TopRow=min(T.TopRow);
BottomRow=max(T.BottomRow);

LeftCol=min(T.LeftCol);
RightCol=max(T.RightCol);

TboxUnion=    array2table([TopRow   BottomRow   LeftCol   RightCol],...
         'VariableNames',{'TopRow','BottomRow','LeftCol','RightCol'});

end