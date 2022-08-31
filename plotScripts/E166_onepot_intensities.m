%% E166_onepot_intensities plot
% script based on E166_smFISH_SpotSize.m
% which is based on E159part1rep2_ampOverRounds_100X_spotSize_plot.m

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
%parentDirForExport='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';

Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E166_smFISH_onepot_100X_Tcond.xlsx'));

plotDir=[parentDir,filesep,'paper/plots/E166_onepot_amplification/'];
plotName='onepot_amplification';

extractedDataDirForAllCond=fullfile(parentDir,filesep,'paper/extractedData/E166_smFISH_onepot_100X/');


% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig4_onepot.xlsx';
outdataSheetSpots='SupFig4_onepotSpots';
outdataSheetCells='SupFig4_onepotCells';
outdataFileForN='N_SupFig4_onepotSpots.xlsx';
outdataSheetForN='N_SupFig4_onepotSpots';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
%% Import Data
NmPerPx=65; % for 100X, 1x1 binning

Tspots=table();
Tcell=table();
uniqueCondIDsCount=0;
for iExperiment=2

    %if iExperiment==1
    %    experimentID='E159part1rep2';
    %    extractedDataDirForAllCond=fullfile(parentDir,filesep,'paper/extractedData/E159part1rep2_ampclick_100X/');
    %elseif iExperiment==2
    experimentID='E166';


    TspotsThisExp=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tspots1All.csv'),'ReadVariableNames', true);
    TspotsThisExp=TspotsThisExp(TspotsThisExp.isGood==1,:); % only spots where isGood==true
    TspotsThisExp.AmpRound(isnan(TspotsThisExp.AmpRound))=0; % make NaN values 0 values
    TspotsThisExp.experimentID=repmat({experimentID},height(TspotsThisExp),1);

    TcellThisExp=array2table(unique([TspotsThisExp.condID,TspotsThisExp.AmpRound,TspotsThisExp.arrayNum,TspotsThisExp.objNum],'rows'),'VariableNames',{'condID','AmpRound','arrayNum','objNum'});
    TcellThisExp.uniqueCellID=[1:height(TcellThisExp)]';
    TcellThisExp.experimentID=repmat({experimentID},height(TcellThisExp),1);
    TcellThisExp=movevars(TcellThisExp,'experimentID','Before',1);


    TspotsThisExp=join(TspotsThisExp,TcellThisExp(:,{'experimentID','condID','arrayNum','objNum','uniqueCellID'}),'Keys',{'experimentID','condID','arrayNum','objNum'});
    TspotsThisExp.sigmaFittedNm=TspotsThisExp.sigmaFitted * NmPerPx;

    TspotsThisExp=movevars(TspotsThisExp,'uniqueCellID','After','AmpRound');

    TspotsThisExp=movevars(TspotsThisExp,'experimentID','Before',1);
    Tspots=[Tspots;TspotsThisExp];
    Tcell=[Tcell;TcellThisExp];
end

clear TspotsThisExp TcellThisExp
% Append unique condID
tempT=unique(Tspots(:,{'experimentID','condID','condName'}));
tempT.uniqueCondID=[1:height(tempT)]';
Tspots=join(Tspots,tempT);
Tcell=join(Tcell,tempT);

assert(height(Tcell)==length(unique(Tspots.uniqueCellID)))

% only clampFISH
TspotsClamp=Tspots(ismember(Tspots.condID,[2 3 4 5]),:);
TcellClamp=Tcell(ismember(Tcell.condID,[2 3 4 5]),:);

TcountsAll=grpstats(TspotsClamp(:,{'uniqueCellID','channel'}),{'uniqueCellID','channel'});
%a=join(Tcell,tempTcounts,'Keys',{'uniqueCellID','channel'})
channelList={'cyb'};
for iChannel=1:length(channelList)
    channel=channelList{iChannel};
    TcountsChannel=TcountsAll(strcmp(TcountsAll.channel,channel),:);
    TcountsChannel=renamevars(TcountsChannel,'GroupCount',['numSpots_',channel]);
    %Tcell2=outerjoin(Tcell,TcountsChannel,'Keys',{'uniqueCellID'});
    TcellClamp=join(TcellClamp,TcountsChannel,'Keys',{'uniqueCellID'});
    %
end

Tspots=TspotsClamp;
Tcell=TcellClamp;
clear TspotsClamp TcellClamp

%% appearance and get data
YLim=[0 500]; % spot count
FigurePosition=[100 100 540 150];
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
%MedianLineColor=[0 0 0];%black
%MedianLineWidth=3;


jitterOn=true; % for cells per cell


numPlotCols=2;

f=figure(1);clf;
f.Position=FigurePosition;
%t=tiledlayout(1,numPlotCols);
t=tiledlayout(1,2);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='none';
%xlabel(t,'Round of amplification','FontSize',8)
%ylabel(t,sprintf('standard deviation\nof gaussian fit of spots (nm)'))
t.YLabel.FontSize=8;
t.XLabel.Color=[0 0 0]; t.YLabel.Color=[0 0 0];


%         tickLabels={'clampFISH 2.0\newline       {\itMITF}\newline    round 1',...
%                     'clampFISH 2.0\newline       {\itMITF}\newline    round 4',...
%                     'clampFISH 2.0\newline       {\itMITF}\newline one-pot #1',...
%                     'clampFISH 2.0\newline       {\itMITF}\newline one-pot #2'};
tickLabels={'  standard\newlineamplification\newline  round 1',...
    '  standard\newlineamplification\newline  round 4',...
    '   one-pot\newlineamplification\newline       #1',...
    '   one-pot\newlineamplification\newline       #2'};

% Tspots
condID_2_dataGroup=[nan 1 2 3 4]';
Tspots.dataGroup=condID_2_dataGroup(Tspots.condID);
numDataGroups=4;

dataGroup_Tspots=Tspots.dataGroup;
amplitudeFitted=Tspots.amplitudeFitted;

channel='cyb';
idxForNorm=all([strcmp(Tspots.experimentID,'E166'),Tspots.condID==2,Tspots.AmpRound==1,strcmp(Tspots.channel,channel)],2);
normalizationAmplitudeFittedValue=median(Tspots.amplitudeFitted(idxForNorm));
assert(length(unique(Tspots.uniqueCondID(idxForNorm)))==1)
amplitudeFittedNorm=Tspots.amplitudeFitted/normalizationAmplitudeFittedValue;

% Tcell

numSpots_Tcell=Tcell.(['numSpots_',channel]);
dataGroup_Tcell=condID_2_dataGroup(Tcell.condID);

% simple data output
TdataOutSpots=Tspots;
TdataOutSpots.amplitudeFittedNorm=amplitudeFittedNorm;
TdataOutSpots.dataGroup=dataGroup_Tspots;
writetable(TdataOutSpots,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheetSpots)
TdataOutCells=Tcell;
TdataOutCells.dataGroup=dataGroup_Tcell;
writetable(TdataOutCells,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheetCells)
% Export Sample Size N for SPOTS
groupingVars={'condID','condName','channel','AmpRound'};
DataVars={'amplitudeFittedNorm'};
ToutNspots=sortrows(grpstats(TdataOutSpots(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'condName','AmpRound'})
writetable(ToutNspots,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet',outdataSheetForN,'WriteMode','overwrite')
%% PLOT #1 SPOT data: spot amplitude (fitted)
axisNum=1;

%    colNum=2;
%    axisNum=numPlotCols*(i-1)+colNum;
ax=nexttile(t,axisNum); hold on; ax.FontSize=8; title(sprintf('spot amplitude'),'FontWeight','normal','FontSize',8);
ax.YLim=[0.3 16];
%ax.YTick=[1 2 4 8 16 32 64 128];
ax.YTick=[0.5 1 2 4 8 16];
vs=violinplot(amplitudeFittedNorm,dataGroup_Tspots,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
arrayViolinColor=repmat({ViolinColor},[numDataGroups,1]);
[vs(:).ViolinColor]=deal(arrayViolinColor{:});
ax.YScale='log';
% delete the whiskers
if removeWhiskers
    for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
end
ylabel(sprintf('normalized amplitude\nof gaussian fit of spots (AU)'),'FontSize',8,'Color',[0 0 0])
% set median circle size
temp=[vs(:).MedianPlot];
[temp.SizeData]=deal(MedianSizeData);
[vs(:).MedianPlot]=deal(temp);
% label median values with text
medianValues=[temp.YData];
medianValuesStr2f=arrayfun(@(x) num2str(x,'%.2f'),medianValues,'UniformOutput', 0);
%medianValuesStr1f=arrayfun(@(x) num2str(x,'%.1f'),medianValues,'UniformOutput', 0);
%medianValuesStr=[medianValuesStr2f(1:3),medianValuesStr1f(4:6)];
medianValuesStr=medianValuesStr2f;
Xlocs=[1:numDataGroups]-0.4; Ylocs=1.3*medianValues;
text(Xlocs,Ylocs,medianValuesStr,'FontSize',6)
ax.XTick=[1:4]; ax.XTickLabel=tickLabels;
ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
%% PLOT #2 CELL data: spots per cell
axisNum=2;
%axisNum=numPlotCols*(i-1)+colNum;
ax=nexttile(t,axisNum); hold on; ax.FontSize=8; title(sprintf('spot count per cell'),'FontWeight','normal','FontSize',8);

% boxplot
%boxplot(numSpots_Tcell,AmpRoundOrder_Tcell,'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
b=boxchart(dataGroup_Tcell,numSpots_Tcell);%'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
b.MarkerStyle='none';
b.BoxFaceColor=BoxFaceColor;
b.WhiskerLineColor=WhiskerLineColor;
%axMedians=findobj(ax,'Type','Line','Tag','Median');
%[axMedians(:).Color]=deal(MedianLineColor);%
%[axMedians(:).LineWidth]=deal(MedianLineWidth);%
%vs=violinplot(numSpots_Tcell,AmpRound_Tcell,'ShowData',true,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
% data

if jitterOn
    rng(0); jitter=0.3*(rand(size(dataGroup_Tcell))-0.5);
    AmpRoundOrder_Tcell_jitter=dataGroup_Tcell+jitter;
else
    AmpRoundOrder_Tcell_jitter=dataGroup_Tcell;
end
sc=scatter(AmpRoundOrder_Tcell_jitter,numSpots_Tcell,'MarkerFaceColor',numSpotsColor,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5,'SizeData',MarkerSize);
ax.XTick=[1:4]; ax.XTickLabel=tickLabels;
ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
ylabel('number of spots per cell','FontSize',8,'Color',[0 0 0])

ylim(ax,YLim);
%set(findobj(t.Children,'Type','Axes'),'YLim',YLim)


exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']), 'ContentType', 'vector')