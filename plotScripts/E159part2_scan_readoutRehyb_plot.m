% E159part2_scan_readoutRehyb_plot
% spot count correlation in cycle 1 vs. cycle 4 (a repeat of cycle 1)



parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E159part2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'paper/experiments/E159part2_scan_Tscan.xlsx'));
extractedDataDir=fullfile(parentDir,'paper/extractedData/E159part2_scan/scan7_withStripAndRepeatData/');
plotDir=fullfile(parentDir,filesep,'paper/plots/E159part2_scan_readoutRehyb/');

cd(parentDir)

% import Tcell (for drug-naive cells)
Tcell=readtable(fullfile(extractedDataDir,'TcellAll2.csv'));
TplotsAll=readtable(fullfile(plotDir,filesep,'E159part2_scan_readoutRehyb_Tplots.xlsx'));

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig14_rehybRep1.xlsx';
outdataSheet='SupFig14_rehybRep1';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
writetable(Tcell,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)


%% figNum=1
figNum=1;
plotOutName='readoutRehybCorrelations';
Tplots=TplotsAll(TplotsAll.figNum==1,:);
plot_regression=true;
plotYequalsX=true;
export_axes_rasters=true;
rasterResolution=900;

markerSize=2;
markerAlpha=0.2;
markerColor=[0.1 0.1 0.5];
f=figure(figNum); clf
figWidth=550;
figHeight=550;
fontSize=8;
set(f,'Position',[900 200 figWidth figHeight],'Units','points')

numArrayRows=max(Tplots.arrayRow);
numArrayCols=max(Tplots.arrayCol);
t=tiledlayout(numArrayRows,numArrayCols);

t.TileSpacing='compact';
t.Padding='compact';

fontSize=8;

arrayPositionStrList=cell(height(Tplots),1);

for iPlot=1:height(Tplots)
    arrayRow=Tplots.arrayRow(iPlot);
    arrayCol=Tplots.arrayCol(iPlot);
    iAx=arrayCol + numArrayCols*rem(arrayRow-1,numArrayRows);
    figure(figNum)
    ax=nexttile(iAx); hold on; axis equal;
    ax.FontSize=fontSize;
    Xfield=Tplots.Xfield{iPlot};
    Yfield=Tplots.Yfield{iPlot};
    Xlabel=Tplots.Xlabel{iPlot};
    Ylabel=Tplots.Ylabel{iPlot};
    
	X=Tcell.(Xfield);
    Y=Tcell.(Yfield);
    
    
    
    % plot
    s=scatter(ax,X,Y,markerSize,markerColor,'filled');
    s.MarkerFaceAlpha=markerAlpha;
    
    % title
    subplotTitle=Tplots.subplotTitle{iPlot};
    subplotTitle=split(subplotTitle,'\n');
    if ~isempty(subplotTitle)
        title(ax,subplotTitle,'FontSize',fontSize)
    end
    
    % xlabel
    Xlabel=Tplots.Xlabel(iPlot);
    xlabel(Xlabel,'Color','k')
    ax.XAxis.Color=[0 0 0];
    % ylabel
    Ylabel=Tplots.Ylabel(iPlot);
    ylabel(Ylabel,'Color','k')
    ax.XAxis.Color=[0 0 0];
    
    % set x,y limit
    gene=Tplots.Gene{iPlot};
    relevantGeneFields=Tcell.Properties.VariableNames(contains(Tcell.Properties.VariableNames,gene));
    maxData=max(Tcell{:,relevantGeneFields},[],'all');
    xlim([-0.05*maxData 1.05*maxData])
    ylim([-0.05*maxData 1.05*maxData])
    if length(ax.XTick)<=length(ax.YTick)
        ax.YTick=ax.XTick;
    else
        ax.XTick=ax.YTick;
    end
    
   

    % linear regression
    mdl=fitlm(X,Y);
    Tcoeff=mdl.Coefficients;
    m=Tcoeff.Estimate('x1'); % slope
    b=Tcoeff.Estimate('(Intercept)'); % intercept
    Rsq=mdl.Rsquared.Ordinary;
    % plot regression
    if plot_regression
        hreg=fplot(ax,@(x) m*x + b);
        hreg.Color='k';
    end

    % plot y=x
    if plotYequalsX
        fplot(ax,@(x) x,'LineStyle',':','Color','k','LineWidth',0.5);
        text(max(xlim)*0.7,max(ylim)*0.8,'y=x','FontSize',fontSize)
    end

    % gene text
    text(max(xlim)*0.05,max(ylim)*0.9,gene,'FontSize',fontSize,'FontAngle','italic')
    

    % plot regression text
    if plot_regression
        %regressionStr=sprintf('R^{2}=%0.3f\ny=%0.2f\\timesx+%.2f',Rsq,m,b);
        regressionStr=sprintf('R^{2}=%0.3f\ny=%0.2fx+%.2f',Rsq,m,b);
        text(max(xlim)*0.65,max(ylim)*0.1,regressionStr,'FontSize',fontSize)
    end
    
    ax.Tag=['ax',num2str(iAx),'_r',num2str(arrayRow),'c',num2str(arrayCol)];
end

exportgraphics(f,fullfile(plotDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector');
%% export everyting but without the datapoints

fNoData=figure; set(fNoData,'Position',get(f,'Position'));
set(fNoData,'Units','points')
copyobj(f.Children,fNoData);
delete(findobj(fNoData.Children.Children,'Type','Scatter'))
delete(findobj(fNoData.Children.Children,'Type','FunctionLine'))
%set(fNoData.Children.Children,'XTick',[],'YTick',[])
%set(fNoData.Children.Children,'XColor','none','YColor','none')

exportgraphics(fNoData,fullfile(plotDir,filesep, [plotOutName,'_NoData.eps']), 'ContentType', 'vector');

%% export rasterized axes (with scatter, lines, ticks but NOT tick labels, axes labels, regression text, or titles)
if export_axes_rasters
    fRasterParent=figure;
    set(fRasterParent,'Position',get(f,'Position'))
    set(fRasterParent,'Units','points')
    copyobj(f.Children,fRasterParent);
    set(fRasterParent.Children.Children,'XTickLabel',[],'YTickLabel',[])
    set(fRasterParent.Children.Children,'XLabel',[],'YLabel',[])
    set(fRasterParent.Children.Children,'Title',[])
    delete(findobj(fRasterParent.Children.Children,'Type','Text'))
end

%output axis rasters to subfolder
if isfolder(fullfile(plotDir,filesep,'axes_raster',filesep))
    rmdir(fullfile(plotDir,filesep,'axes_raster',filesep),'s')
end
if ~isfolder(fullfile(plotDir,filesep,'axes_raster',filesep))
    mkdir(fullfile(plotDir,filesep,'axes_raster',filesep))
end

for iAx=1:length(fRasterParent.Children.Children)
        axRaster=fRasterParent.Children.Children(iAx);
        axPositionStr=axRaster.Tag;
        exportgraphics(axRaster,fullfile(plotDir,filesep,'axes_raster',filesep,[plotOutName,'_',axPositionStr,'.jpg']), 'ContentType', 'image','Resolution',rasterResolution);
end
    