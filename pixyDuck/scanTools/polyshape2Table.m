function T=polyshape2Table(tileIDs,polyvec)

assert(size(tileIDs,2)==1) % is column vector
assert(size(tileIDs,1)==length(polyvec));


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

T=array2table([tileIDs, rStart,rEnd,cStart,cEnd],'VariableNames',{'tileID','TopRow','BottomRow','LeftCol','RightCol'});


end