function [rStart,rEnd,cStart,cEnd]=polyshape2RowColStartEnd(polyvec)
% second dimension of vertices interpreted as row
% first dimension of vertices interpreted as column

numPoly=length(polyvec);

rStart=nan(numPoly,1);
rEnd=nan(numPoly,1);
cStart=nan(numPoly,1);
cEnd=nan(numPoly,1);

for i=1:numPoly
   
    vertices=polyvec(i).Vertices;
    
    rStart(i)=min(vertices(:,2));
    rEnd(i)=max(vertices(:,2));
    
    cStart(i)=min(vertices(:,1));
    cEnd(i)=max(vertices(:,1));
    
end

end