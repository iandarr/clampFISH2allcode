% E157part1_strip_plot
% readout stripping

%% Import data

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E157part1rep2_strip_Tcond.csv'));

cd(parentDir)

% Get paired (before/after strip) Tcell data
extractedDataDirForAllCond='paper/extractedData/E157part1rep2_strip/';
TcellAll=readtable(fullfile(parentDir,filesep,extractedDataDirForAllCond,filesep,'TcellAll2.csv'));

plotOutDir=fullfile(parentDir,filesep,'paper/plots/E157part1rep2_strip/');
%% read Table of plot specifications and generate X,G for plotting
cd(parentDir)
Tplots=readtable(fullfile(plotOutDir,filesep,'Tplots.xlsx'));

plotGroupIDs=Tplots.plotGroupID';

Ydata=[]; Xgroup=[];
YdataPaired=[];XgroupPaired=[];
for i=1:length(plotGroupIDs)
    plotGroupID=plotGroupIDs(i);
    
    x_read=(2*i-1);
    x_strip=(2*i);
    
    scanID=Tplots.scanID(plotGroupID);
    numericChannelPrefix=Tplots.numericChannel{plotGroupID};
    
    indTcellPlotPaired=TcellAll.scanID==scanID;
    numPoints=sum(indTcellPlotPaired);
%     if ~height(unique(TcellAll(indTcellPlotPaired,{[numericChannelPrefix,'_gene']}),'rows'))==1
%         error('this condition has multiple gene values') 
%     end
%     assert(all(TcellAll.s_condID(indTcellPlotPaired)==condID_strip))
    
    Ydata_read=TcellAll.(['read_',numericChannelPrefix,'_numSpots'])(indTcellPlotPaired);
    Ydata_strip=TcellAll.(['strip_',numericChannelPrefix,'_numSpots'])(indTcellPlotPaired);
    Ydata=[Ydata;Ydata_read;Ydata_strip];
    Xgroup=[Xgroup;repmat(x_read,numPoints,1);repmat(x_strip,numPoints,1)];
    
    YdataPaired=[YdataPaired;[Ydata_read,Ydata_strip]];
    XgroupPaired=[XgroupPaired;[repmat(x_read,numPoints,1),repmat(x_strip,numPoints,1)]];
end

%% boxplots of counts before and after strip
close all
figureWidthPoints=600;
figureHeightPoints=300;


for iPlot=1:2
    if iPlot==1
        f=figure(1); 
        plotType='full';
        set(f,'Position',[1500 600 figureWidthPoints figureHeightPoints],'Units','points')
        plotOutName='ReadoutStrip_full';
    else
        f=figure(2);
        plotType='zoom';
        set(f,'Position',[1500 200 figureWidthPoints figureHeightPoints],'Units','points')
        plotOutName='ReadoutStrip_zoom';
    end
    clf;

ax=gca;hold on;
readColor=[0.5 0.8 0.6];
stripColor=[0.9 0.5 0.2];

numGroupsRead=length(unique(XgroupPaired(:,1)));
positionsVect_centers=1:numGroupsRead;
positionsVect_read=positionsVect_centers-0.2;
positionsVect_strip=positionsVect_centers+0.2;
set(ax,'XLim',[positionsVect_centers(1)-1,positionsVect_centers(end)+1])

outlierMarkerSize=2;
boxplot(ax,YdataPaired(:,1),XgroupPaired(:,1),'positions',positionsVect_read,'Symbol','o','Colors',readColor,'BoxStyle','filled','Width',1,'FactorGap',0);
hOutliers=findobj(ax.Children(1),'Tag','Outliers');
set(hOutliers,'MarkerEdgeColor','none','MarkerFaceColor',readColor,'MarkerSize',outlierMarkerSize)

YLim=ax.YLim;
boxplot(ax,YdataPaired(:,2),XgroupPaired(:,2),'positions',positionsVect_strip,'Symbol','o','Colors',stripColor,'BoxStyle','filled','Width',1,'FactorGap',0);
hOutliers=findobj(ax.Children(1),'Tag','Outliers');
set(hOutliers,'MarkerEdgeColor','none','MarkerFaceColor',stripColor,'MarkerSize',outlierMarkerSize)


set(ax,'YLim',YLim)
if strcmp(plotType,'zoom')
    set(ax,'YLim',[-10 500])
end
set(ax,'XLim',[positionsVect_centers(1)-1,positionsVect_centers(end)+1])
ax.XTick=positionsVect_centers;
ax.XTickLabels=positionsVect_centers;

XTickLabels=Tplots.gene;
%XTickLabels(9:10)=join([XTickLabels(9:10),repmat({'{\n}X'},2,1)],'');
ax.XTickLabels=XTickLabels;

% box width
boxWidth=7;
hBoxes=findobj(ax,'Type','Line','Tag','Box');
set(hBoxes,'LineWidth',boxWidth)

% median width
medianWidth=0.2
hMedians=findobj(gca,'Tag','Median');
set(hMedians,'Color','k','LineWidth',0.5)
for i=1:length(hMedians)
    XData=hMedians(i).XData;
    XData=mean(XData)+[-medianWidth/2, medianWidth/2]; 
    hMedians(i).XData=XData;
end

fontSize=8;
set(ax,'FontSize',fontSize)
ylabel('Spots per cell')

% whiskers
hWhiskers=findobj(gca,'Tag','Whisker');


%ax.YLim=[-10 200]
text(0.02*ax.XLim(2), 0.95*ax.YLim(2),'Before Strip','Color',readColor,'FontSize',fontSize)
text(0.02*ax.XLim(2), 0.9*ax.YLim(2),'After Strip','Color',stripColor,'FontSize',fontSize)

exportgraphics(f, fullfile(plotOutDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector'); % with spots
end



%%  before-after scatterplot for a given gene
% %numericChannel='A594';
% gene='FN1';
% 
% readCondID=Tplots.condID_readout(strcmp(Tplots.gene,gene));
% stripCondID=Tplots.condID_strip(strcmp(Tplots.gene,gene));
% numericChannel=Tplots.numericChannel{strcmp(Tplots.gene,gene)};
% 
% %indTcellPlotPaired=strcmp(TcellPlotPaired.([numericChannel,'_gene']),gene);
% indTcellPlotPaired=TcellPlotPaired.r_condID==readCondID;
% readData=TcellPlotPaired{indTcellPlotPaired,(['r_',numericChannel,'_numSpots'])};
% stripData=TcellPlotPaired{indTcellPlotPaired,(['s_',numericChannel,'_numSpots'])};
% X=readData;
% Y=stripData;
% 
% figure(2)
% ph=plot(X,Y,'.');
% ph.MarkerSize=14;
% 
% ax=gca;
% 
% % title
% title(sprintf('%s',gene))
% 
% % LABELS
% labelsCellStr=replace(join([strtrim(cellstr(num2str(TcellPlotPaired.arrayNum(indTcellPlotPaired)))),...
%                                             strtrim(cellstr(num2str(TcellPlotPaired.objNum(indTcellPlotPaired))))],2),' ','-');
% cx=0.01*ax.XLim(2);
% cy=0.01*ax.YLim(2);
% text(X+cx,Y+cy,labelsCellStr)
% 
% xlabel('Spots BEFORE Strip')
% ylabel('Spots AFTER Strip')
% ax.FontSize=16
% 
