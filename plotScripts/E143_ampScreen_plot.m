%% E143_ampScreen_plot
% amplifier screen


%% Go to parent folder
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)

%% Read in data
Tcond=readtable('paper/experiments/E143_amplifierScreen_Tcond.xlsx');

extractedDataDir='paper/extractedData/E143_screen/60X';

condIDList=[16:30]';
%condIDList=[1:45]';

Tspots=readtable(fullfile(extractedDataDir,filesep,'TspotsAll.csv'));
TspotsClampTopN=readtable(fullfile(extractedDataDir,filesep,'TspotsClampTopN.csv'));
Tcell=readtable(fullfile(extractedDataDir,filesep,'TcellAll3_smfishN_clampAvgIntNabove20.csv'));
TcondResults=readtable(fullfile(extractedDataDir,filesep,'TcondResults.csv'));

%% boxplots of EGFR and GFP spot intensities
plotOutDir='paper/plots/E143_amplifier_screen';
plotOutName='amplifierScreen_spotIntensityBoxplots';

plotList={'GFP','EGFR'};

TspotsClampTopN2=join(TspotsClampTopN,Tcond(:,{'condID','primary','amplifierSeries'}),'Keys','condID');

% for simple plot data output - ALL data points (not just medians)
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig5_ampscreen.xlsx';
outdataSheet='SupFig5_ampscreen';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
writetable(TspotsClampTopN2,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)
% Export Sample Size N for spots
outdataFileForN='N_SupFig5_ampscreen.xlsx';
groupingVars={'condID','channel','primary','amplifierSeries'};
DataVars={'intensities'};
ToutNspots=sortrows(grpstats(TspotsClampTopN2(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'condID','amplifierSeries'});
writetable(ToutNspots,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','N_SupFig5_ampscreen','WriteMode','overwrite')


for iPlot=1:2
    geneName=plotList{iPlot};
f=figure(iPlot); clf;

figureWidthPoints=300;
figureHeightPoints=300;
set(f,'Position',[1500 200 figureWidthPoints figureHeightPoints],'Units','points')



showOutliers=false;

axesLabelFontSize=8;
%condIDInPlotOrder=sort(unique(TspotsClampTopN.condID),'ascend');

if strcmp(geneName,'GFP')
    condIDInPlotOrder=[16:30]';
elseif strcmp(geneName,'EGFR')
    condIDInPlotOrder=[31:45]';
end

if showOutliers
    outlierSymbol='k.';
else
    outlierSymbol='';
end

TspotsClampTopNToPlot=TspotsClampTopN2(ismember(TspotsClampTopN2.condID,condIDInPlotOrder),:);

hb=boxplot(TspotsClampTopNToPlot.intensities,TspotsClampTopNToPlot.amplifierSeries,'OutlierSize',4,'Symbol',outlierSymbol,'Widths',0.6,'PlotStyle','compact','MedianStyle','line');%,'Notch','on');%,'PlotStyle','compact','Jitter',0.2)
ax=gca;
boxwhiskerColor=[182 219 255]/255; % light blue
set(hb(1,:),"Color",boxwhiskerColor) % whisker color
set(hb(2,:),"Color",boxwhiskerColor) % box color
set(hb(3,:),"Color",'k') % median lines black

% number of outliers
numTotalSpots=height(TspotsClampTopNToPlot);
axOutliers=findobj(ax,'Type','Line','Tag','Outliers');
numOutliers=length([axOutliers.YData]);
fprintf('plot %i (%s): %i outliers out of %i points total (%.3f %%) \n',iPlot,geneName,numOutliers,numTotalSpots,100*numOutliers/numTotalSpots)


ax.FontSize=axesLabelFontSize;
ax.YAxis.Exponent=0;


% whiskers
axWhiskers=findobj(ax,'Type','Line','Tag','Whisker');
set(axWhiskers,'LineWidth',1) % whisker width

% set YLim
if ~showOutliers
    set(ax,'YLim',[0 1.05*max([axWhiskers.YData])])
end

% box width
boxWidth=6;
axBoxes=findobj(ax,'Type','Line','Tag','Box');
set(axBoxes,'LineWidth',boxWidth)

% median line length
medianLength=0.45;
axMedians=findobj(ax,'Type','Line','Tag','Median');
for i=1:length(axMedians)
    axMedians(i).XData=[-medianLength/2, medianLength/2] + mean(axMedians(i).XData);
end

if strcmp(geneName,'GFP')
    title('GFP clampFISH 2.0 spot intensities')
else
    title('{\itEGFR} clampFISH 2.0 spot intensities')
end




ylabel('spot intensity (AU)')
%set(ax,'XTickLabelRotation',0)
xlabel('amplifier probe set','FontSize',axesLabelFontSize)

% export figure
exportgraphics(f, fullfile(plotOutDir,filesep, [plotOutName,'_',geneName,'.eps']), 'ContentType', 'vector'); % with spots
end

%% EGFR vs. GFP
plotOutDir='paper/plots/E143_amplifier_screen';
plotOutName='amplifierScreen_EGFR_vs_GFP';

% for simple plot data output - MEDIANS only
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig6_ampmedians.xlsx';
outdataSheet='SupFig6_ampmedians';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.


f=figure(3); clf
figureWidthPoints=380;
figureHeightPoints=275;
set(f,'Position',[1500 200 figureWidthPoints figureHeightPoints],'Units','points')

plot_regression_line=true;

GFP_condID=16:30;
EGFR_condID=31:45;
TcondX=TcondResults(GFP_condID,:);
TcondY=TcondResults(EGFR_condID,:);
assert(all([TcondX.amplifierSeries==TcondY.amplifierSeries],1))
assert(all(strcmp(TcondResults.primary(GFP_condID),'GFP')) && all(strcmp(TcondResults.primary(EGFR_condID),'EGFR')))

plotMetric='median';
markerSize=14;
labelFontSize=8;
regressionFontSize=8;
axesLabelFontSize=8;

X=TcondX.(plotMetric);
Y=TcondY.(plotMetric);

% for simple plot data output - MEDIANS only
Tdataout=array2table([[1:15]',X,Y],'VariableNames',{'amplifier set','GFP primaries','EGFR primaries'});
writetable(Tdataout,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)


plot(X,Y,'.','MarkerSize',markerSize); hold on
ax=gca;
ax.XLim=[0 9000];
ax.YLim=[0 9000];
ax.FontSize=axesLabelFontSize;
axis equal;

if ax.YLim(2)>ax.XLim(2),ax.XLim=ax.YLim;end

title(sprintf("%s spot intensity of each amplifier probe set",plotMetric))

lax=legend();
    if plot_regression_line
        [r,m,b]=regression(X',Y'); % use data before the jitter
        if abs(m)<1
            Xreg=ax.XLim';
            Yreg=m*Xreg+b;
        else
            Yreg=ax.YLim';
            Xreg=(Yreg-b)/m;
        end            
        plot(Xreg,Yreg,'-','Color','k')
        lax.String=lax.String(1:end-1);
        text(ax.XLim(2)*.60,Yreg(2)*.6,sprintf('R^2=%.3f\ny=%.2fx+%.0f',r^2,m,b),'FontSize',regressionFontSize)
    end
legend off;
xlabel([sprintf('%s spot intensity (AU) with 10 ',plotMetric),'GFP primary probes'])
ylabel([sprintf('%s spot intensity (AU) with 30 ',plotMetric),'{\itEGFR} primary probes'])

XlabelOffset=-450*(Y>=m*X+b) + 50*(Y<m*X+b);
YlabelOffset=150   *(Y>=m*X+b) + -150*(Y<m*X+b);
text(X+XlabelOffset,Y+YlabelOffset,cellstr(num2str(TcondX.amplifierSeries)),'FontSize',labelFontSize)

ax.XTick=ax.YTick;
exportgraphics(f, fullfile(plotOutDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector'); % with spots
