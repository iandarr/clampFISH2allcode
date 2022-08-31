function TinnerBox=getInnerBoxFromPolyvec(tilesPolyvec)

poly=getInnerRectPoly2(tilesPolyvec);

XYvertices=[poly.Vertices];
% x is column, y is row

assert(size(XYvertices,1)==4)

TopRow=min(XYvertices(:,2));
assert(sum(XYvertices(:,2)==TopRow)==2)

BottomRow=max(XYvertices(:,2));
assert(sum(XYvertices(:,2)==BottomRow)==2)

LeftCol=min(XYvertices(:,1));
assert(sum(XYvertices(:,1)==LeftCol)==2)

RightCol=max(XYvertices(:,1));
assert(sum(XYvertices(:,1)==RightCol)==2)


TinnerBox=array2table([ TopRow   BottomRow   LeftCol   RightCol],...
     'VariableNames', {'TopRow','BottomRow','LeftCol','RightCol'});




end