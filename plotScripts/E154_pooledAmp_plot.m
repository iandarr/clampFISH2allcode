% E154_pooledAmp_plot
% pooled amplification
%% Import Tcell and Tspots from all folders
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)

Tcond=readtable('paper/experiments/E154_pooled_amp_Tcond.xlsx');

extractedDataDir='paper/extractedData/E154_pooled_amp/s6 60X p2';

Tspots=readtable(fullfile(extractedDataDir,filesep,'TspotsClampTopK3_1.csv'));
%Tspots=readtable(fullfile(extractedDataDir,'TspotsClampTopK_3.csv')); % with fitted

Tspots2=join(Tspots,Tcond(:,{'condID','amplifierSeries','amplifiersIncluded'}),'Keys','condID');
%Tspots2.amplifierSeries=categorical(Tspots2.amplifierSeries,[1,3,5,6,7,9,10,12,14,15]);

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig7_pooledAmp.xlsx';
outdataSheet='SupFig7_pooledAmp';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
writetable(Tspots2,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)
% Export Sample Size N for SPOTS
groupingVars={'condID','channel','amplifierSeries','amplifiersIncluded'};
DataVars={'intensities'};
ToutNspots=sortrows(grpstats(Tspots2(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'amplifierSeries','amplifiersIncluded'})
outdataFileForN='N_SupFig7_pooledAmp.xlsx';
writetable(ToutNspots,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','N_SupFig7_pooledAmp','WriteMode','overwrite')




%% conditions to plot
condToPlot=25:44;
showOutliers=false;
fontSize=8;
plotOutDir='paper/plots/E154_pooled_amp';
plotOutName='AmplifierIntensities_PooledVsAlone';

if showOutliers
    outlierSymbol='k.';
else
    outlierSymbol='';
end

f=figure(1); clf
figureWidthPoints=540;
figureHeightPoints=300;
set(f,'Position',[1500 200 figureWidthPoints figureHeightPoints],'Units','points')
ax=gca;
ax.FontSize=fontSize;

aloneColor=[0.5 0.8 0.6];
pooledColor=[0.9 0.5 0.2];

numAmpSeries=length(unique(Tspots2.amplifierSeries));

ampSeriesList=[1,3,5,6,7,9,10,12,14,15];
centerXlist=1:10;
ampSeries2centerX=containers.Map(ampSeriesList,centerXlist);


posAlone=centerXlist-0.15;
posPool=centerXlist+0.15;

set(ax,'XLim',[centerXlist(1)-1,centerXlist(end)+1])

idxAlone=strcmp(Tspots2.amplifiersIncluded,'alone');

YAlone=Tspots2.intensities(idxAlone);
ampSetAlone=Tspots2.amplifierSeries(idxAlone);
boxplot(ax,YAlone,ampSetAlone,'positions',posAlone,'Symbol',outlierSymbol,'Colors',aloneColor,'BoxStyle','filled','Width',1,'FactorGap',0);
hOutliersAlone=findobj(ax.Children(1).Children(1),'Tag','Outliers');
numOutliersAlone=length([hOutliersAlone(:).XData]);
fprintf('ALONE data:  %i outliers out of %i points total\n',numOutliersAlone,length(ampSetAlone))

hold on;

idxPool=strcmp(Tspots2.amplifiersIncluded,'pool');
YPool=Tspots2.intensities(idxPool);
ampSetPool=Tspots2.amplifierSeries(idxPool);
boxplot(ax,YPool,ampSetPool,'positions',posPool,'Symbol',outlierSymbol,'Colors',pooledColor,'BoxStyle','filled','Width',1,'FactorGap',0);

hOutliersPool=findobj(ax.Children(1).Children(1),'Tag','Outliers');
numOutliersPool=length([hOutliersPool(:).XData]);
fprintf('POOLED data: %i outliers out of %i points total\n',numOutliersPool,length(ampSetPool))

% box width
boxWidth=6;
axBoxes=findobj(ax,'Type','Line','Tag','Box');
set(axBoxes,'LineWidth',boxWidth)

% median line length
medianLength=0.18;
axMedians=findobj(ax,'Type','Line','Tag','Median');
for i=1:length(axMedians)
    axMedians(i).XData=[-medianLength/2, medianLength/2] + mean(axMedians(i).XData);
end
% median line color and width
[axMedians(:).Color]=deal([1 1 1]);%

% whiskers
axWhiskers=findobj(ax,'Type','Line','Tag','Whisker');
set(axWhiskers,'LineWidth',1) % whisker width

% set YLim
if ~showOutliers
    set(ax,'YLim',[0 1.05*max([axWhiskers.YData])])
end

%set(ax,'YLim',YLim)
%set(ax,'YLim',[YLim(1) 1000])
set(ax,'XLim',[centerXlist(1)-1,centerXlist(end)+1])
ax.XTick=centerXlist;
ax.YAxis.Exponent=0;

title(ax,'spot intensities for each amplifier set when used alone vs. in a pool of 10 amplifier sets')
ylabel(ax,sprintf('GFP clampFISH 2.0 spot intensities (AU)\nof top 10,000 spots in 40 cells'))
xlabel(ax,'amplifier set')

% legend

text(0.02*ax.XLim(2), 0.95*ax.YLim(2),'alone','Color',aloneColor,'FontSize',fontSize,'FontWeight','bold')
text(0.02*ax.XLim(2), 0.9*ax.YLim(2),'pooled','Color',pooledColor,'FontSize',fontSize,'FontWeight','bold')

exportgraphics(f, fullfile(plotOutDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector'); % with spots
