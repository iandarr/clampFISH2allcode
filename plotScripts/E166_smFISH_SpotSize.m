%% E166_smFISH_SpotSize plot
% script based on E159part1rep2_ampOverRounds_100X_spotSize_plot
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
%parentDirForExport='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';

Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E166_smFISH_onepot_100X_Tcond.xlsx'));
Tplots=readtable(fullfile(parentDir,filesep,'paper/plots/E159part1rep2_E166_spotSize/Tplots_E166.xlsx'));
Tplots.condIDs=cellfun(@(x) str2num(x),Tplots.condIDs,'UniformOutput',false);
Tplots.experimentIDs=cellfun(@(x) replace(split(x,',')',' ',''),Tplots.experimentIDs,'UniformOutput',false);
Tplots.channels=cellfun(@(x) replace(split(x,',')',' ',''),Tplots.channels,'UniformOutput',false);

plotDir=[parentDir,filesep,'paper/plots/E159part1rep2_E166_spotSize/'];
plotName='smFISH_vs_clamp2_SpotSize100X';

dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SourceData_ExtDataFig6_SpotSizeVsSmFISH.xlsx';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
outdataFileForN='N_ExtDataFig6_SpotSizeVsSmFISH.xlsx';
%% Import Data
NmPerPx=65; % for 100X, 1x1 binning

Tspots=table();
Tcell=table();
uniqueCondIDsCount=0;
for iExperiment=1:2
    
    if iExperiment==1
        experimentID='E159part1rep2';
        extractedDataDirForAllCond=fullfile(parentDir,filesep,'paper/extractedData/E159part1rep2_ampclick_100X/');
    elseif iExperiment==2
        experimentID='E166';
        extractedDataDirForAllCond=fullfile(parentDir,filesep,'paper/extractedData/E166_smFISH_onepot_100X/');
    end

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
tempT=unique(Tspots(:,{'experimentID','condID'}));
tempT.uniqueCondID=[1:height(tempT)]';
Tspots=join(Tspots,tempT);
Tcell=join(Tcell,tempT);

% % count spots per dye and put in Tcell
% for i=1:height(Tplots)
%     dye=Tplots.dyeToPlot{i};
%     condIDs=Tplots.condIDs(i);
%     for iCond=length(condIDs)
%         
%         condID=condIDs(iCond);
% 
%         TtempCell=grpstats(Tspots(strcmp(Tspots.channel,channel),{'uniqueCellID','sigmaFitted'}),"uniqueCellID");
%         TtempCell=renamevars(TtempCell,'GroupCount',['numSpots_',channel]);
%         Tcell=join(Tcell,TtempCell(:,{'uniqueCellID',['numSpots_',channel]}));
%     end
% end
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

%% appearance and get data
colLabelNames={'Atto 488 labeled';...
    'Atto 647N labeled'};
YLim=[0 250];
FigurePosition=[50 50 540 120];
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


numPlotCols=max(Tplots.plotCol);

f=figure(1);clf;
f.Position=FigurePosition;
%t=tiledlayout(1,numPlotCols);
t=tiledlayout(1,8);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='none';
%xlabel(t,'Round of amplification','FontSize',8)
ylabel(t,sprintf('standard deviation\nof gaussian fit of spots (nm)'))
t.YLabel.FontSize=8;
t.XLabel.Color=[0 0 0]; t.YLabel.Color=[0 0 0];

for i=1:numPlotCols
    
    if i==1
        axisNum=1;
        axisSpan=[1 3];
        tickLabels={'smFISH\newline {\itUBC}','clampFISH 2.0\newline      {\itUBC}\newline    round 1','clampFISH 2.0\newline      {\itUBC}\newline    round 4'};
    else
        axisNum=4;
        axisSpan=[1 5];
        tickLabels={'smFISH\newline{\itTOP2A}','clampFISH 2.0\newline       {\itMITF}\newline    round 1','clampFISH 2.0\newline       {\itMITF}\newline    round 4','clampFISH 2.0\newline       {\itMITF}\newline    round 1','clampFISH 2.0\newline       {\itMITF}\newline    round 4'};
    end
    labelName=colLabelNames{i};

    rowTplots=find(Tplots.plotCol==i);
    condIDs=Tplots.condIDs{rowTplots};
    experimentIDs=Tplots.experimentIDs{rowTplots};
    channels=Tplots.channels{rowTplots};
    
    numDataGroups=length(condIDs);
    TspotsSubset=table();
    for iDataGroup=1:numDataGroups
        condID=condIDs(iDataGroup);
        experimentID=experimentIDs{iDataGroup};
        channel=channels{iDataGroup};
        idxSubsetTspots=all([strcmp(Tspots.experimentID,experimentID),Tspots.condID==condID,strcmp(Tspots.channel,channel)],2); % subset
        assert(length(unique(Tspots.uniqueCondID(idxSubsetTspots)))==1) % doubel check all spots come from one unique condition
        
        TspotsSubsetThisDataGroup=Tspots(idxSubsetTspots,:);
        TspotsSubsetThisDataGroup.dataGroup=repmat(iDataGroup,height(TspotsSubsetThisDataGroup),1);
        TspotsSubset=[TspotsSubset;TspotsSubsetThisDataGroup];
    end

    % Export simple data used for plot
    groupingVars={'experimentID','condID','condName','AmpRound','channel'};
    DataVars='sigmaFittedNm';
    ToutN=sortrows(grpstats(TspotsSubset(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'experimentID','condName'});
    
    if i==1
        writetable(TspotsSubset,fullfile(dirForDataExport,filesep,outdataFile),'Sheet','UBC_Atto488')
        
        % export N spots
        writetable(ToutN,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','UBC_Atto488','WriteMode','overwritesheet');

    elseif i==2
        writetable(TspotsSubset,fullfile(dirForDataExport,filesep,outdataFile),'Sheet','MITF_Atto647N')
        
        % export N spots
        writetable(ToutN,fullfile(dirForDataExport,filesep,outdataFileForN),'Sheet','MITF_Atto647N','WriteMode','overwritesheet');
    end

    % Tspots
    dataGroup=TspotsSubset.dataGroup;
    sigmaFittedNm=TspotsSubset.sigmaFittedNm;

%     if i==1
%         idxForNorm=all([strcmp(Tspots.experimentID,'E159part1rep2'),Tspots.AmpRound==1,strcmp(Tspots.channel,'gfpb')],2);
%     else
%         idxForNorm=all([strcmp(Tspots.experimentID,'E159part1rep2'),Tspots.AmpRound==1,strcmp(Tspots.channel,'cyb')],2);
%     end
%     normalizationAmplitudeFittedValue=median(Tspots.amplitudeFitted(idxForNorm));
%     assert(length(unique(Tspots.uniqueCondID(idxForNorm)))==1)
%     amplitudeFittedNorm=TspotsSubset.amplitudeFitted/normalizationAmplitudeFittedValue;

    % Tcell
    %numSpots_Tcell=Tcell.(['numSpots_',channel]);
    %AmpRound_Tcell=Tcell.AmpRound;

    % SPOT data: spot size
    
    %axisNum=numPlotCols*(i-1)+colNum;
    %axisNum=colNum;
    ax=nexttile(t,axisNum,axisSpan); hold on; ax.FontSize=8; title(sprintf('%s spot size',labelName));
    vs=violinplot(sigmaFittedNm,dataGroup,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    arrayViolinColor=repmat({ViolinColor},[numDataGroups,1]);
    [vs(:).ViolinColor]=deal(arrayViolinColor{:});
    % set median circle size
    temp=[vs(:).MedianPlot];
    [temp.SizeData]=deal(MedianSizeData);
    [vs(:).MedianPlot]=deal(temp);
    % label median values with text
    medianValues=[temp.YData];
    medianValuesStr=arrayfun(@(x) num2str(x,'%.0f'),medianValues,'UniformOutput', 0);
    Xlocs=[1:numDataGroups]-0.3; Ylocs=medianValues +15;
    text(Xlocs,Ylocs,medianValuesStr,'FontSize',6)
    
    ax.XAxis.TickLabels=tickLabels;
    ax.XAxis.Color=[0 0 0]; ax.YAxis.Color=[0 0 0];
    %ax.XAxis.TickLabelColor=[0 0 0]; ax.YAxis.TickLabelColor=[0 0 0];
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
    end
    %ylabel(sprintf('standard deviation\nof gaussian fit of spots (nm)'))

%     % SPOT data: spot amplitude (fitted)
%     colNum=2;
%     axisNum=numPlotCols*(i-1)+colNum;
%     ax=nexttile(t,axisNum); hold on; ax.FontSize=8; title(sprintf('%s',labelName)); %title(sprintf('%s amplitude of gaussian fit',labelName));
%     ax.YLim=[0.4 200];
%     ax.YTick=[1 2 4 8 16 32 64 128];
%     vs=violinplot(amplitudeFittedNorm,dataGroup,'ShowData',ShowData,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
%     arrayViolinColor=repmat({ViolinColor},[numDataGroups,1]);
%     [vs(:).ViolinColor]=deal(arrayViolinColor{:});
%     ax.YScale='log';
%     % delete the whiskers
%     if removeWhiskers
%         for iViolin=1:length(vs), vs(iViolin).WhiskerPlot.LineStyle='none'; end
%     end
%     ylabel(sprintf('normalized amplitude\nof gaussian fit of spots (AU)'))
%     % set median circle size
%     temp=[vs(:).MedianPlot];
%     [temp.SizeData]=deal(MedianSizeData);
%     [vs(:).MedianPlot]=deal(temp);
%     % label median values with text
%     medianValues=[temp.YData];
%     medianValuesStr2f=arrayfun(@(x) num2str(x,'%.2f'),medianValues,'UniformOutput', 0);
%     %medianValuesStr1f=arrayfun(@(x) num2str(x,'%.1f'),medianValues,'UniformOutput', 0);
%     %medianValuesStr=[medianValuesStr2f(1:3),medianValuesStr1f(4:6)];
%     medianValuesStr=medianValuesStr2f;
%     Xlocs=[1:numDataGroups]; Ylocs=1.2*medianValues;
%     text(Xlocs,Ylocs,medianValuesStr,'FontSize',6)
% 
%     % CELL data: spots per cell
%     colNum=3;
%     axisNum=numPlotCols*(i-1)+colNum;
%     ax=nexttile(t,axisNum); hold on; ax.FontSize=8; %title(sprintf('%s mplitude of gaussian fit',labelName));
%     convertVect=[1 2 nan 3 nan 4 nan 5 nan 6]';
%     AmpRoundOrder_Tcell=convertVect(AmpRound_Tcell);
%     % boxplot
%     %boxplot(numSpots_Tcell,AmpRoundOrder_Tcell,'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
%     b=boxchart(AmpRoundOrder_Tcell,numSpots_Tcell);%'Symbol','','BoxStyle','filled','Colors',boxplotColor,'Width',1,'MedianStyle','line') %'Width',1,'FactorGap',0 % outliers already visible due to the scatter
%     b.MarkerStyle='none';
%     b.BoxFaceColor=BoxFaceColor;
%     b.WhiskerLineColor=WhiskerLineColor;
%     %axMedians=findobj(ax,'Type','Line','Tag','Median');
%     %[axMedians(:).Color]=deal(MedianLineColor);%
%     %[axMedians(:).LineWidth]=deal(MedianLineWidth);%
%     %vs=violinplot(numSpots_Tcell,AmpRound_Tcell,'ShowData',true,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
%     % data
% 
%     if jitterOn
%         rng(0); jitter=0.3*(rand(size(AmpRound_Tcell))-0.5);
%         AmpRoundOrder_Tcell_jitter=AmpRoundOrder_Tcell+jitter;
%     else
%         AmpRoundOrder_Tcell_jitter=AmpRoundOrder_Tcell;
%     end
%     sc=scatter(AmpRoundOrder_Tcell_jitter,numSpots_Tcell,'MarkerFaceColor',numSpotsColor,'MarkerEdgeColor','none','MarkerFaceAlpha',0.5,'SizeData',MarkerSize);
%     ax.XTick=[1:6];ax.XTickLabel={'1','2','4','6','8','10'};
% 
%     ylabel('number of spots per cell')
end
set(findobj(t.Children,'Type','Axes'),'YLim',YLim)


exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']), 'ContentType', 'vector')