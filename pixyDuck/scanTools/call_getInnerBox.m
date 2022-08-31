% test_getInnerBox

load(fullfile('/Volumes/IAND_05/E159part2_scan/scan/scan1_WellA1',filesep,'ScanObject.mat'))

thisRound=1;
res=1;
tilesPolyvec=Scan.Rounds(thisRound).fullTiles.rscan.res(res).polyvec;

%%
getInnerRectPoly2(tilesPolyvec,true)
%%
load('/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/polyvecProbIn.mat');
polyvec=polyIn;
getInnerRectPoly2(polyvec,true)
%%
fh=figure(1); clf;
ax=axes(fh); 
poly1 = polyshape([0 0 1 1],[1 0 0 1]);
poly2 = polyshape([0.75 1.25 1.25 0.75],[0.25 0.25 0.75 0.75]);
poly3 = polyshape([1.1 1.1 1.3 1.3],[.6 1.5 1.5 0.6]);
polyvec=[poly1;poly2; poly3];
%plot(polyvec); hold on;
%%
poly=union(polyvec);
plot(poly); hold on
%%
innerBoxPoly=getInnerRectPoly(polyvec,true)


%%
fh=figure(1); clf;
ax=axes(fh); 
%tileRange=1:10;
pvec1=plot(ax,tilesPolyvec);
hold on; set(ax,'YDir','reverse')
[pvec1.FaceColor]=deal([0 1 1]); 
labels=arrayfun(@num2str,[1:length(tilesPolyvec)]','UniformOutput',0);
[X,Y]=tilesPolyvec.centroid;
X=X'; Y=Y';
text(X,Y,labels,'FontSize',8)
ax.XLim=fh.Children(1).XLim; ax.YLim=fh.Children(1).YLim;
% %pvec2=plot(ax,Scan.Rounds(2).polyvec.regscan.FullImgs(tileRange));
% axis(ax,'equal');
% 
% xlim =
% 
%    1.0e+04 *
% 
%    -0.0361    7.2203
% 
% 
% ylim =
% 
%    1.0e+03 *
% 
%     3.4340    3.9622

%%
innerBoxPoly=getInnerBoxFromPolyvec4(tilesPolyvec,4)
%%
poly=RowColStartEnd2polyshape(Tbox.TopRow,Tbox.BottomRow,Tbox.LeftCol,Tbox.RightCol);

pvec2=plot(poly)
[pvec2.FaceColor]=deal([1 0 0]);

%%
figure(1); clf; ax=gca;
plot(Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).polyvec)
hold on;axis(ax,'equal');set(ax,'YDir','reverse');
plot(Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res(iRes).polyvec)

%%
figure(1); clf; ax=gca; set(ax,'YDir','reverse'); hold on;
for thisRound=1:10%[1 5]
    polyvec=Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res.polyvec;
    plot(polyvec)
    pause(1)
end
%%

plot(Scan.regions.fullTiles.innerBoxIntersect.rscan.res.polyvec)