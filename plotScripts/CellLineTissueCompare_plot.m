% CellLineTissueCompare_plot.m
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
extractedDataDir=fullfile(parentDir,'paper/extractedData/CellLineTissueComparison',filesep);

Tspots=readtable(fullfile(extractedDataDir,filesep,'Tspots.csv'));
plotDir=[parentDir,filesep,'paper/plots/CellLineTissueComparison/'];
plotName='CellLineVsTissueSpotAmplitude';

grpstats(Tspots(:,{'amplitude','intensity','sigma','sampleType','BackgLevel'}),'sampleType',{'mean','median','min'})

% for simple plot data output - ALL data points (not just medians)
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig18_CellLineTissueCompare.xlsx';
outdataSheet='SupFig18_CellLineTissueCompare';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
% Export Sample Size N for SPOTS
outdataFileForN='N_SupFig18_CellLineTissueCompare.xlsx';

%%
% figure(1)
% X=Tspots.intensity;
% Y=Tspots.amplitude;
% g=Tspots.sampleType;
% gscatter(X,Y,g);
% xlabel('intensity'); ylabel('amplitude')
%ylim([0 0.1])

% %%
% figure(2)
% X=Tspots.sigma;
% Y=Tspots.amplitude;
% g=Tspots.sampleType;
% gscatter(X,Y,g);
% xlabel('sigma');
% ylabel('amplitude')
TspotsSorted = sortrows(Tspots,{'sampleType'},'ascend');
sampleTypes=unique(TspotsSorted.sampleType,'stable');

YLim=[0 7];
FigurePosition=[50 50 400 300];
% appearance of violins
ViolinAlpha=1;
ViolinColor=[182 219 255]/255; % light blue
%numSpotsColor=[36 255 36]/255;

% appearanch of violin plots
removeWhiskers=true; % they don't really add anything.
ShowData=false; % show all points

% appearance of median marker
MedianMarker='o'; %default 'o'
MedianLineWidth=0.02; %default 0.5
MedianSizeData=14; % default 36
MedianColor=[1 1 1]; %[1 1 1]=white

% appearance of interquartile range box
BoxWidth=0.02;
BoxColor=[0.4 0.4 0.4];

% spots per cell: scatter + boxplot on
%boxplotColor=[0 105 208]/255;
MarkerSize=10;
BoxFaceColor=[0 105 208]/255;
BoxFaceColor=[0.5 0.5 0.5];
WhiskerLineColor=[0.5 0.5 0.5];
numSpotsColor=[0 110 218]/255;


numPlotCols=1;
%sampleTypes={'CellLine','TissueFreshFrozen','TissueFFPE'};

tickLabels=replace(sampleTypes,{'CellLine','TissueFreshFrozen','TissueFFPE'},...
    {'      Cell Line\newline(WM989 A6-G3)', '   Tissue, Fresh Frozen\newline(WM989-A6-G3-Cas9-5a3)','    Tissue, FFPE\newline(CPLX WM4505-1)'})

f=figure(1);clf;
f.Position=FigurePosition;
%t=tiledlayout(1,numPlotCols);
t=tiledlayout(1,1);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='none';
%xlabel(t,'Round of amplification','FontSize',8)
ylabel(t,sprintf('amplitude of gaussian fit of spots, normalized'),'FontSize',8)
t.YLabel.FontSize=8;
t.XLabel.Color=[0 0 0]; t.YLabel.Color=[0 0 0];

normAmplitude=median(TspotsSorted.amplitude(strcmp(TspotsSorted.sampleType,'CellLine')));

% for simple plot data output - ALL data points (not just medians)
Tdataout=[TspotsSorted,array2table(TspotsSorted.amplitude/normAmplitude,'VariableNames',{'amplitudeNorm'})];
writetable(Tdataout,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)
% Export Sample Size N for SPOTS
groupingVars={'sampleType','channel'};
DataVars={'amplitudeNorm'};
ToutNspots=sortrows(grpstats(Tdataout(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'sampleType'})

writetable(ToutNspots,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','N_SupFig18_CellTissue','WriteMode','overwrite')


for i=1:numPlotCols
    axisNum=i;
    axisSpan=[1 1];
    % Tspots
    dataGroup=TspotsSorted.sampleType;
    numDataGroups=length(unique(dataGroup));
    amplitude=TspotsSorted.amplitude;
    amplitudeNorm=amplitude/normAmplitude;
    % SPOT data: spot size
    ax=nexttile(t,axisNum,axisSpan); hold on; ax.FontSize=8;
    %title(sprintf('spot amplitude'));
    vs=violinplot(amplitudeNorm,dataGroup,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    arrayViolinColor=repmat({ViolinColor},[numDataGroups,1]);
    [vs(:).ViolinColor]=deal(arrayViolinColor{:});
    % set median circle size
    temp=[vs(:).MedianPlot];
    [temp.SizeData]=deal(MedianSizeData);
    [vs(:).MedianPlot]=deal(temp);
    % label median values with text
    medianValues=[temp.YData];
    medianValuesStr=arrayfun(@(x) num2str(x,'%.2f'),medianValues,'UniformOutput', 0);
    Xlocs=[1:numDataGroups]-0.17; Ylocs=medianValues;
    text(Xlocs,Ylocs,medianValuesStr,'FontSize',8)

    ax.XAxis.TickLabels=tickLabels;
    ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
    %ax.XAxis.TickLabelColor=[0 0 0]; ax.YAxis.TickLabelColor=[0 0 0];
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
    end
end
    set(findobj(t.Children,'Type','Axes'),'YLim',YLim)


exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']), 'ContentType', 'vector')
