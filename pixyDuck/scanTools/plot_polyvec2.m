%%
clear Scan
scanDir='/Volumes/IAND_07/rawData/E157part1rep2_stripping/registeredstitches/SetA_resistant';
load(fullfile(scanDir,filesep,'ScanObject.mat'))
%%
clf;
fh=figure(1);
ax=axes(fh); 

tileRange=1:10;

thisRound=1;
res=1;
labels=arrayfun(@num2str,[1:length(polyvec)]','UniformOutput',0);
X=Scan.Rounds(1).fullTiles.rscan.res.T.CenterCol;
Y=Scan.Rounds(1).fullTiles.rscan.res.T.CenterRow;
text(ax,X,Y,labels,'FontSize',8)

hold on; set(ax,'YDir','reverse')
polyvec=Scan.Rounds(thisRound).fullTiles.rscan.res(res).polyvec;
pvec1=plot(ax,polyvec,'FaceAlpha',0.1);
[pvec1.FaceColor]=deal([1 0 0]); 
%hold on
%pvec2=plot(ax,Scan.Rounds(2).polyvec.regscan.FullImgs(tileRange));
axis(ax,'equal');
% plot round 2
% clf;
% fh=figure(1);
% ax=axes(fh);
thisRound=2;
res=1;
polyvec=Scan.Rounds(thisRound).fullTiles.rscan.res(res).polyvec;
pvec2=plot(ax,polyvec,'FaceAlpha',0.1);
[pvec2.FaceColor]=deal([0 0 1]); 
hold on;

%% plot boxes
%plot(Scan.regions.fullTiles.innerBoxIntersect.rscan.res.polyvec)
%plot(Scan.regions.fullTiles.outerBoxUnion.rscan.res.polyvec,'FaceAlpha',0)
plot(Scan.Rounds(1).regions.fullTiles.innerBox.rscan.res.polyvec,'FaceAlpha',0.2)

%%
fh2=figure
ax=axes(fh2); 
hold on; set(ax,'YDir','reverse')

colorList=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
%for thisRound=1:4
innerPoly=Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res(res).polyvec;
outerPoly=Scan.Rounds(thisRound).regions.fullTiles.outerBox.rscan.res(res).polyvec;

innerIntersectPoly=Scan.regions.fullTiles.innerBoxIntersect.rscan.res(res).polyvec;

%plot(outerPoly,'FaceColor','k','FaceAlpha',0.1)
hIn(thisRound)=plot(innerPoly,'EdgeColor',colorList(thisRound,:),'FaceColor',colorList(thisRound,:),'FaceAlpha',0);
%end
innerBoxIntersect=Scan.regions.fullTiles.innerBoxIntersect.rscan.res(1).polyvec;
outerBoxIntersect=Scan.regions.fullTiles.outerBoxIntersect.rscan.res(1).polyvec;
plot(innerBoxIntersect,'LineStyle','none','EdgeColor','k','FaceColor','k','FaceAlpha',0.1);
plot(outerBoxIntersect,'LineStyle','-','EdgeColor','k','FaceColor','k','FaceAlpha',0);
%%
hInInter=plot(innerIntersectPoly,'FaceColor','k','FaceAlpha',0.3)
%%
X=Scan.Rounds(1).fullTiles.rstage.T.centerX;
Y=Scan.Rounds(1).fullTiles.rstage.T.centerY;
hold on
text(X,Y,labels,'FontSize',8)
%hold on
%pvec2=plot(ax,Scan.Rounds(2).polyvec.regscan.FullImgs(tileRange));
axis(ax,'equal');

