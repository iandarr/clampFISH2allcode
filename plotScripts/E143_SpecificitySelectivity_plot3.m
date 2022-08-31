% E143_SpecificitySelectivity_plot
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E143_amplifierScreen_Tcond.xlsx'));

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
extractedDataDir=fullfile(parentDir,filesep,'paper/extractedData/E143_SensitivitySelectivity/');
% Import Data
Tres=readtable(fullfile(extractedDataDir,filesep,['TsensitivitySelectivity.csv']));

plotDir=fullfile(parentDir,filesep,'paper/plots/E143_SensitivitySelectivity/');

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig2_coloc.xlsx';
outdataSheet='SupFig2_coloc';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
writetable(Tres(Tres.isShifted==false,:),fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)
%% global options
% Switch plotType between FullData and SimpleData 
%plotType='FullData';
plotType='SimpleData';

separateFigureForLegendOn=true;
figNumForLegend=2;
plotAnyError=false;
asteriskToMarkKeptAmplifierSeries=true;
keptAmplifierSeries=[1 3 5 6 7 9 10 12 14 15];
plotLegend=false;
iPlotForLegend=3;
% lines to plot
minDenom=100;
%condIDstoPlot=[(1:15)',(16:30)'];
condIDstoPlot=[16:30]'; % with primaries; without primaries will be taken automatically
%condIDstoPlot=[31:45]'; % with primaries; without primaries will be taken automatically
%condIDstoPlot=[16:45]'; % with primaries; without primaries will be taken automatically
%condIDstoPlot=[16:19]'; % with primaries; without primaries will be taken automatically
Tlines=table();
darkBlue=[0,0,139]/255;
lightBlue=[182 219 255]/255;
darkOrange=[255, 140, 0]/255;
lightOrange=[255, 197, 127]/255;

Tlines.colorSpec=[darkBlue;darkBlue;lightBlue;darkOrange;darkOrange;lightOrange];
Tlines.lineStyle=repmat({'-','-','-'}',2,1);
Tlines.plotError=repmat([true false false]',2,1);
Tlines.isShifted=repmat([false true false]',2,1);
Tlines.withPrimary=repmat([true true false]',2,1);
Tlines.legendLabels=repmat({'(+) primary probes','(+) primary probes, 10px shifted','(-) primary probes'}',2,1);
Tlines.dataToPlot={'clamp','clamp','clamp','smFISH','smFISH','smFISH'}';
Tlines.lineWidth=repmat([3,1,1]',2,1);
Tlines.markerType=repmat({'o','^','v'}',2,1);
Tlines.markerSize=repmat([4 3 3]',2,1);
if strcmp(plotType,'FullData')
    % Tlines is full height
elseif strcmp(plotType,'SimpleData')
    Tlines=Tlines([1,4],:);
else
    error('plotType mst be FullData or SimpleData')
end
FigurePosition=[200 200 560 600];

MarkerSize=10;

figNum=1;
f=figure(figNum);clf;
f.Position=FigurePosition;
numPlots=length(condIDstoPlot);
numPlotCols=3;
numPlotRows=ceil(numPlots/numPlotCols);
t=tiledlayout(numPlotRows,numPlotCols);
if separateFigureForLegendOn
    fLegend=figure(figNumForLegend);
    fLegend.Position=[200 20 560 100]; hold off;
end
%t=tiledlayout(1,1);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='compact';
xlabel(t,'intensity threshold (fluorescence units)','FontSize',8)
ylabel(t,'fraction of spots colocalized','FontSize',8)


for iPlot=1:numPlots


    ax=nexttile(t,iPlot); hold on;
    ax.XScale='log';
    for iLine=1:height(Tlines)
        isShifted=Tlines.isShifted(iLine);
        withPrimary=Tlines.withPrimary(iLine);
        dataToPlot=Tlines.dataToPlot{iLine};
        plotError=Tlines.plotError(iLine);

        % formatting
        lineStyle=Tlines.lineStyle{iLine};
        colorSpec=Tlines.colorSpec(iLine,:);
        lineWidth=Tlines.lineWidth(iLine);
        markerType=Tlines.markerType{iLine};
        markerSize=Tlines.markerSize(iLine);
        if withPrimary
            condID=condIDstoPlot(iPlot);
            amplifierSeries=unique(Tres.amplifierSeries(Tres.condID==condID));
        else
            condIDofWithPrimary=condIDstoPlot(iPlot);
            amplifierSeries=unique(Tres.amplifierSeries(Tres.condID==condIDofWithPrimary));
            condID=unique(Tres.condID(all([Tres.amplifierSeries==amplifierSeries,strcmp(Tres.primary,'NoPrimaries')],2)));
        end
        idxAll=all([Tres.isShifted==isShifted,Tres.condID==condID],2); % all denominator (total spots) values
        if strcmp(dataToPlot,'smFISH')
            idx=idxAll & (Tres.numSmFISHTot>=minDenom);
            X=Tres.threshold(idx);
            Ynum=Tres.numSmFISHColoc(idx);
            Ydenom=Tres.numSmFISHTot(idx);
        elseif strcmp(dataToPlot,'clamp')
            idx=idxAll & (Tres.numClampTot>=minDenom);
            X=Tres.threshold(idx);
            Ynum=Tres.numClampColoc(idx);
            Ydenom=Tres.numClampTot(idx);
        else
            error('unkown case')
        end

        if strcmp(dataToPlot,'clamp') && withPrimary
            XdispRange=[min(X),max(X)];
            primaryGene=unique(Tres.primary(idx));assert(length(primaryGene)==1);primaryGene=primaryGene{1};
            if asteriskToMarkKeptAmplifierSeries  && ismember(amplifierSeries,keptAmplifierSeries)
                asteriskStr='*';
            else
                asteriskStr='';
            end
            %titleStr=sprintf('%s, amplifier set %i %s',primaryGene,amplifierSeries,asteriskStr);
            titleStr=sprintf('amplifier set %i %s',amplifierSeries,asteriskStr);
            %Xticks=X(linspace(1,length(X),min(4,length(X))));
            Xticks=X;
        end
        [Ymid,Yrange]=binofit(Ynum,Ydenom);
        %Y=Ynum./Ydenom;
        plot(X,Ymid,markerType,'LineStyle',lineStyle,'Color',colorSpec,'LineWidth',lineWidth,'MarkerSize',markerSize);
        
        if (iPlot==1) & separateFigureForLegendOn
            figure(figNumForLegend);% figure for legend only
            plot(X,Ymid,markerType,'LineStyle',lineStyle,'Color',colorSpec,'LineWidth',lineWidth,'MarkerSize',markerSize);
            hold on;
            figure(figNum); % back to normal figure
        end
    
        if plotAnyError
            if plotError
                %plot(X,Yrange(:,1),'-');
                %plot(X,Yrange(:,2),'-');
                X2=[X;flipud(X)];
                Y2=[Yrange(:,1);flipud(Yrange(:,2))];
                fill(X2,Y2,colorSpec,'FaceAlpha',0.3,'EdgeColor','none')
            end
        end
        %legendLabel=Tlines.legendLabels{iLine};

    end % end iLines loop
    title(ax,titleStr,'FontWeight','normal','FontSize',8)
    ax.FontSize=8;
    ax.XLim=XdispRange;
    ax.YLim=[0 1];
    ax.XTick=Xticks;
    YTicks=[0 0.25 0.50 0.75 1]';
    ax.YTick=YTicks;
    ax.YTickLabel=cellstr(num2str(YTicks,'%0.2f'));
    if ismember(iPlot,iPlotForLegend) & plotLegend
        lh=legend(ax,Tlines.legendLabels);
    end
    % legend plot in second figure
end
if separateFigureForLegendOn
    figure(figNumForLegend);
    axLegendPlot=gca;
    [axLegendPlot.Children(:).Visible]=deal('off');
    axLegendPlot.Box='off';
    axLegendPlot.XAxis.Visible='off';
    axLegendPlot.YAxis.Visible='off';
    %legend on;
    lhand=legend(Tlines.legendLabels);
    axLegendPlot.XLim=[0 1];axLegendPlot.YLim=[0 1];
    lhand.Position=[0.1 0.1 0.8 0.8];
    figure(figNum)
end
exportgraphics(f,fullfile(plotDir,filesep,['E143_SensitivitySpecificity_',plotType,'.eps']))
exportgraphics(fLegend,fullfile(plotDir,filesep,['E143_SensitivitySpecificity_',plotType,'Legend.eps']))

