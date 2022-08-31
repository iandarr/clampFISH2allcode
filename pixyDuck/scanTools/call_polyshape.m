
% tiles1
% rStart=[2 4]';
% rEnd=  [9 20]';
% cStart=[1 10]';
% cEnd=  [4 5]';

% last doesnt overlap
rStart=[2 4 23]';
rEnd=  [9 20 36]';
cStart=[1 10 8]';
cEnd=  [4 5 16]';

polyVect_tiles=RowColStartEnd2polyshape(rStart,rEnd,cStart,cEnd);

% rStart=[1]';
% rEnd=  [20]';
% cStart=[3]';
% cEnd=  [30]';
% polyVect_subregion=RowColStartEnd2polyshape(rStart,rEnd,cStart,cEnd);
figure(1)
plot([polyVect_tiles,polyVect_subregion])
set(gca,'YDir','reverse')
XLim=get(gca,'XLim');
YLim=get(gca,'YLim');
%
[rStartOut,rEndOut,cStartOut,cEndOut]=polyshape2RowColStartEnd(polyVect_tiles);
rs=[rStart,rStartOut]
re=[rEnd,rEndOut]
cs=[cStart,cStartOut]
ce=[cEnd,cEndOut]
%%

% call intersect
polyvec_intersect=intersect(polyVect_tiles,polyVect_subregion)
figure(2)
plot(polyvec_intersect)
set(gca,'YDir','reverse')
set(gca,'XLim',XLim)
set(gca,'YLim',YLim)

[rStart,rEnd,cStart,cEnd]=polyshape2RowColStartEnd(polyvec_intersect);


%% overlaps
overlaps(polyVect_tiles,polyVect_subregion)