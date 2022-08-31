function polyoutVectOut=RowColStartEnd2polyshape(rStart,rEnd,cStart,cEnd)
% rStart=2;
% rEnd=  9;
% cStart=1;
% cEnd=4;
% 
% rStart=[2 4]';
% rEnd=  [9 20]';
% cStart=[1 10]';
% cEnd= [4 5]';

numPoly=length(rStart);

%polyoutCellArray=cell(numPoly,1);
polyoutVectOut=polyshape.empty([0,numPoly]);
for i=1:numPoly
    x=[cStart(i) cEnd(i) cEnd(i) cStart(i)];
    y=[rStart(i) rStart(i) rEnd(i) rEnd(i)];
    pgon=polyshape(x,y);
    %polyoutCellArray{i}=pgon;
    polyoutVectOut(i)=pgon;
    %plot(pgon1)
    %set(gca,'YDir','reverse')
end

%polyoutVectOut=[polyoutCellArray{:}]';

end