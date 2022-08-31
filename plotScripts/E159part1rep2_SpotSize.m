%% E159part1rep2_SpotSize
% was called E159part1rep2_ampOverRounds_100X_spotSize_plot
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
%parentDirForExport='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';

Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E159part1rep2_ampclick_Tcond.xlsx'));

plotDir=[parentDir,filesep,'paper/plots/E159part1rep2_E166_spotSize/'];
plotName='clamp2_OverRounds_SpotSize100X';

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SourceData_ExtDataFig5_SpotSizeWithAmp.xlsx';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
%% theoretical spot sizes
%
% Gaussian approximations of fluorescence microscope point-spread function models
% Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin
% https://opg-optica-org.proxy.library.upenn.edu/ao/fulltext.cfm?uri=ao-46-10-1819&id=130945#t7


lambdaEmission=[535 667] % midpoint of emission filter band
%lambdaEmission=[535+15 667+15] % upper wavelength of emission filter band
% paraxial
NA=1.45;
theoreticalGaussianFitStd_Paraxial=0.21*lambdaEmission/NA

% nonparaxial
n=1.518; % Refractive index. Nikon Type N immersion oil MXA22165
alpha=asin(NA/n);
k=(2*pi()./lambdaEmission); % wavenumber [1/nm]

theoreticalGaussianFitStd_NonParaxial=...
    (1./(n*k))*...
    ((4 - 7*cos(alpha)^(3/2) + 3*cos(alpha)^(7/2))/...
     (7*(1 - cos(alpha)^(3/2)))       )^(-1/2)

%% appearance
channelsToPlot={'gfpb','cyb'};
labelNames={'{\itUBC}, amplifier set 9, readout probe in Atto 488',...
    '{\itMITF}, amplifier set 12, readout probe in Atto 647N'};

FigurePosition=[50 50 540 275];
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

NmPerPx=65;

% Import Data
extractedDataDirForAllCond=fullfile(parentDir,filesep,'paper/extractedData/E159part1rep2_ampclick_100X/');

Tspots=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tspots1All.csv'),'ReadVariableNames', true);
%Tspots=Tspots(Tspots.isGood==1,:); % only spots where isGood==true
Tcell=array2table(unique([Tspots.condID,Tspots.AmpRound,Tspots.arrayNum,Tspots.objNum],'rows'),'VariableNames',{'condID','AmpRound','arrayNum','objNum'});
Tcell.uniqueCellID=[1:height(Tcell)]';
Tspots=join(Tspots,Tcell(:,{'condID','arrayNum','objNum','uniqueCellID'}),'Keys',{'condID','arrayNum','objNum'});
Tspots.sigmaFittedNm=Tspots.sigmaFitted * NmPerPx;

Tspots=movevars(Tspots,'uniqueCellID','After','AmpRound');


%TtempCell
%Tcell=join(Tcell,)

jitterOn=true; % for cells per cell



for i=1:length(channelsToPlot)
    channel=channelsToPlot{i};
    TtempCell=grpstats(Tspots(strcmp(Tspots.channel,channel),{'uniqueCellID','sigmaFitted'}),"uniqueCellID");
    TtempCell=renamevars(TtempCell,'GroupCount',['numSpots_',channel]);
    Tcell=join(Tcell,TtempCell(:,{'uniqueCellID',['numSpots_',channel]}));
end

% Export simple data used for plot
TdataSpotsOut=Tspots(ismember(Tspots.channel,channelsToPlot),:);
TdataSpotsOut.GeneAmpsetDye=replace(TdataSpotsOut.channel,{'gfpb','cyb'},{'UBCset9Atto488','MITFset12Atto647N'});
TdataSpotsOut=movevars(TdataSpotsOut,'GeneAmpsetDye','After','channel');
writetable(TdataSpotsOut,fullfile(dirForDataExport,filesep,outdataFile),'Sheet','spots')
writetable(Tcell,fullfile(dirForDataExport,filesep,outdataFile),'Sheet','cells')
% Export Sample Size N for SPOTS
groupingVars={'condID','condName','GeneAmpsetDye'};
DataVars={'amplitudeFitted'};
ToutNspots=sortrows(grpstats(TdataSpotsOut(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'GeneAmpsetDye'})
outdataFileForN='N_ExtDataFig5_SpotSizeWithAmp.xlsx';
writetable(ToutNspots,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','N_ExtDataFig5_spots','WriteMode','overwrite')

% Export Sample Size N for CELLS
groupingVars={'condID','AmpRound'};
DataVars={'uniqueCellID'};
ToutNcells=sortrows(grpstats(Tcell(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'AmpRound'})
ToutNcells.median_uniqueCellID=[];
writetable(ToutNcells,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','N_ExtDataFig5_cells','WriteMode','overwrite')

numPlotRows=length(channelsToPlot);
numConditions=length(unique(Tcell.condID));
f=figure(1);clf;
f.Position=FigurePosition;
numPlotCols=3;
t=tiledlayout(2,numPlotCols);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='none';
xlabel(t,'Round of amplification','FontSize',8)
t.XLabel.Color=[0 0 0];t.YLabel.Color=[0 0 0];

for i=1:numPlotRows
    channel=channelsToPlot{i};
    labelName=labelNames{i};
    TspotsSubset=Tspots(strcmp(Tspots.channel,channel),:); % subset by channel

    % Tspots
    AmpRound=TspotsSubset.AmpRound;
    sigmaFittedNm=TspotsSubset.sigmaFittedNm;
    amplitudeFittedNorm=TspotsSubset.amplitudeFitted/median(TspotsSubset.amplitudeFitted(TspotsSubset.AmpRound==1));

    % Tcell
    numSpots_Tcell=Tcell.(['numSpots_',channel]);
    AmpRound_Tcell=Tcell.AmpRound;

    %% SPOT data: spot size
    colNum=1;
    axisNum=numPlotCols*(i-1)+colNum;
    ax=nexttile(t,axisNum); hold on; ax.FontSize=8; %title(sprintf('%s standard deviation of gaussian fit',labelName));
    ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
    %ax.XAxis.TickLabelColor=[0 0 0]; ax.YAxis.TickLabelColor=[0 0 0];
    vs=violinplot(sigmaFittedNm,AmpRound,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    arrayViolinColor=repmat({ViolinColor},[numConditions,1]);
    [vs(:).ViolinColor]=deal(arrayViolinColor{:});
    % set median circle size
    temp=[vs(:).MedianPlot];
    [temp.SizeData]=deal(MedianSizeData);
    [vs(:).MedianPlot]=deal(temp);
    % label median values with text
    medianValues=[temp.YData];
    medianValuesStr=arrayfun(@(x) num2str(x,'%.0f'),medianValues,'UniformOutput', 0);
    Xlocs=[1:6] + -0.6; Ylocs=medianValues +15;
    text(Xlocs,Ylocs,medianValuesStr,'FontSize',6)

    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
    end
    ylabel(sprintf('standard deviation\nof gaussian fit of spots (nm)'))

    %% SPOT data: spot amplitude (fitted)
    colNum=2;
    axisNum=numPlotCols*(i-1)+colNum;
    ax=nexttile(t,axisNum); hold on; ax.FontSize=8; title(sprintf('%s',labelName)); %title(sprintf('%s amplitude of gaussian fit',labelName));
    ax.YLim=[0.4 200];
    ax.YTick=[1 2 4 8 16 32 64 128];
    ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
    %ax.XAxis.TickLabelColor=[0 0 0]; ax.YAxis.TickLabelColor=[0 0 0];
    vs=violinplot(amplitudeFittedNorm,AmpRound,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    arrayViolinColor=repmat({ViolinColor},[numConditions,1]);
    [vs(:).ViolinColor]=deal(arrayViolinColor{:});
    ax.YScale='log';
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
    end
    ylabel(sprintf('normalized amplitude\nof gaussian fit of spots (AU)'))
    % set median circle size
    temp=[vs(:).MedianPlot];
    [temp.SizeData]=deal(MedianSizeData);
    [vs(:).MedianPlot]=deal(temp);
    % label median values with text
    medianValues=[temp.YData];
    medianValuesStr2f=arrayfun(@(x) num2str(x,'%.2f'),medianValues,'UniformOutput', 0);
    medianValuesStr1f=arrayfun(@(x) num2str(x,'%.1f'),medianValues,'UniformOutput', 0);
    medianValuesStr=[medianValuesStr2f(1:3),medianValuesStr1f(4:6)];
    Xlocs=[1:6]+ -0.85; Ylocs=1.2*medianValues;
    text(Xlocs,Ylocs,medianValuesStr,'FontSize',6)

    %% CELL data: spots per cell
    colNum=3;
    axisNum=numPlotCols*(i-1)+colNum;
    ax=nexttile(t,axisNum); hold on; ax.FontSize=8; %title(sprintf('%s mplitude of gaussian fit',labelName));
    ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
    convertVect=[1 2 nan 3 nan 4 nan 5 nan 6]';
    AmpRoundOrder_Tcell=convertVect(AmpRound_Tcell);
    % boxplot
    %boxplot(numSpots_Tcell,AmpRoundOrder_Tcell,'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
    b=boxchart(AmpRoundOrder_Tcell,numSpots_Tcell);%'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
    b.MarkerStyle='none';
    b.BoxFaceColor=BoxFaceColor;
    b.WhiskerLineColor=WhiskerLineColor;
    %axMedians=findobj(ax,'Type','Line','Tag','Median');
    %[axMedians(:).Color]=deal(MedianLineColor);%
    %[axMedians(:).LineWidth]=deal(MedianLineWidth);%
    %vs=violinplot(numSpots_Tcell,AmpRound_Tcell,'ShowData',true,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    % data

    if jitterOn
        rng(0); jitter=0.3*(rand(size(AmpRound_Tcell))-0.5);
        AmpRoundOrder_Tcell_jitter=AmpRoundOrder_Tcell+jitter;
    else
        AmpRoundOrder_Tcell_jitter=AmpRoundOrder_Tcell;
    end
    sc=scatter(AmpRoundOrder_Tcell_jitter,numSpots_Tcell,'MarkerFaceColor',numSpotsColor,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5,'SizeData',MarkerSize);
    ax.XTick=[1:6];ax.XTickLabel={'1','2','4','6','8','10'};

    ylabel('number of spots per cell')

end


exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']), 'ContentType', 'vector')