% E159part2rep2_wStrip_wInt_OpticalCrowdingPlot
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tscan.xlsx'));

extractedDataDir=fullfile(parentDir,'paper/extractedData/E159part2rep2_scan/scan3_withStripAndRepeatData_withMeans',filesep);
scanIDsList=3;

Tchannels=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scanWithStrip_ChannelDesc.xlsx'));
Tchannels=Tchannels(~ismember(Tchannels.channelName,{'Merged_Rcyto_YFP_DAPI','Merged_Rcyto_YFP_DAPI_cp_masks','R123_maxYFP'}),:);
TchannelsForPlot1=Tchannels(1==Tchannels.forCrowdingPlot1,:);
TchannelsForPlot2=Tchannels(1==Tchannels.forCrowdingPlot2,:);
Tcell=readtable(fullfile(extractedDataDir,'TcellAll2_withMeans.csv')); % only has UBC QC checks for 3 of the imaging rounds (R1,R2,R3)

plotDir=fullfile(parentDir,filesep,'paper/plots/E159part2_OpticalCrowding/');

% readout dye and exposure time in tbale
TchannelInfo=table();
TchannelInfo.channelName=TchannelsForPlot1.channelName;
temp=split(TchannelInfo.channelName,'_');
TchannelInfo.exposureTime=temp(:,3);
TchannelInfo.gene=temp(:,4);
channel2Dye={'YFP','Atto 488';'CY3','Cy3';'A594','Alexa Fluor 594';'CY5','Atto 647N'};
[~,ind]=ismember(temp(:,2),channel2Dye(:,1));
TchannelInfo.readoutDye=channel2Dye(ind,2);

gene2amplifierSet={'UBC','set 9';'WNT5A','set 7';'DDX58','set 6'; 'AXL','set 3'; 'NGFR','set 1';'FN1','set 5';'EGFR','set 15';'ITGA3','set 10';'MMP1','set 14';'MITF','set 12'};
[~,ind]=ismember(temp(:,4),gene2amplifierSet(:,1));
TchannelInfo.amplifierSet=gene2amplifierSet(ind,2);
%
FigurePosition=[200 200 640 700];

plotSparse=false; % whether to downsample by sparsityFactor. Speeds up plotting.
sparsityFactor=20; % downsample points if plotSparse==true by this factor

MarkerSize=2;
plotWithDensity=true;
markerColorUniform=[0,0,139]/255; % dark blue
markerAlpha=1;

uniformMaxY=true;
maxY=0.52;
%yPercentile=99.95; % use to determine upper Y lim
%maxYPercentileFactor=1.7; % multiply by point defined by yPercentile to get upper Y lim
uniformMaxX=false;
maxX=5000; %for uniformMaxX


plotWithRegressionOnAllData=false;

textYspacing=0.08; %proportion of y range. Only used when textInAxis==true

textInAxis=false;

densityColorLabel='normalized density of datapoints';

%% Plot 1: Spot count vs. signal
TchannelsForThisPlot=TchannelsForPlot1;
f=figure(1);clf;
f.Position=FigurePosition;
numPlots=height(TchannelsForThisPlot);
numPlotCols=3;
numPlotRows=ceil(numPlots/numPlotCols);
t=tiledlayout(numPlotRows,numPlotCols);
%t=tiledlayout(1,1);
t.TileIndexing='columnmajor';
t.TileSpacing='compact';
t.Padding='normal';
xlabel(t,'mean background-subtracted fluorescence (AU)','FontSize',8)
ylabel(t,'spot count per \mum^{2}','FontSize',8)


fprintf('------------ plot 1 ------------\n')
rng(0);

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig11_CountVsSig.xlsx';
outdataSheet='SupFig11_CountVsSig';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
TdataoutPlot=Tcell(:,1:9);

for iPlot=1:numPlots
    %ax=nexttile(t,iPlot); % this weirdly has a bug when specifying the tile with columnmajor indexing, see https://www.mathworks.com/matlabcentral/discussions/highlights/134368-new-in-r2021a-improvements-to-tiled-chart-layout/40538
    ax=nexttile(t);
    ax.FontSize=8; ax.XLabel.FontSize=8; ax.YLabel.FontSize=8; hold on;
    fprintf('plotting iPlot=%i\n',iPlot)
    varNameForSpots=TchannelsForThisPlot.channelName{iPlot};
    varNameForMeans=['sig_',varNameForSpots];

    X=Tcell.(varNameForMeans);
    Y=Tcell.(varNameForSpots)./Tcell.Area_Um2;
%     Xmean=mean(X); Xmedian=median(X); Xskew=Xmean/Xmedian;
%     Ymean=mean(Y); Ymedian=median(Y); Yskew=Ymean/Ymedian;
%     SkewRatio=Yskew/Xskew;
%     fprintf('iPlot=%2i %25s  Xmean=%10.2f  Xmed=%10.2f  Ymean=%5.4f  Ymed=%5.4f  Xskew=%8.2f  Yskew=%8.2f  SkewRatio=%10.2f\n',iPlot,varNameForSpots,Xmean,Xmedian,Ymean,Ymedian,Xskew,Yskew,SkewRatio)
    

if plotSparse
        % downsample the points by sparsityFactor
        rng(0);
        indSparse=sort(randperm(length(X),round(length(X)/sparsityFactor)))';
        X=X(indSparse);
        Y=Y(indSparse);
    end
    LM=fitlm(X,Y);

    TdataoutPlot.([varNameForSpots,'_mean'])=X;
    TdataoutPlot.([varNameForSpots,'_spot'])=Y;

    if plotWithDensity
        %Vq=ksdensity([X,Y],[X,Y]); % SLOW - use pts to evaluate kernel density are the same locations as the points themselves
        % faster way to use ksdensity is to use meshgrid of points and later interpolate density at the x,y points
        gridx=unique(prctile(X,0:1:100));
        gridy=unique(prctile(Y,0:1:100));
        [Xi,Yi]=meshgrid(gridx,gridy);
        [fdens,~]=ksdensity([X,Y],[Xi(:),Yi(:)]);
        Vq=interp2(Xi,Yi,reshape(fdens,size(Xi)),X,Y);%interpolate at the points of X,Y
        markerColor=Vq/max(Vq); % 1 is highest density
    else
        markerColor=markerColorUniform;
    end
    %scatter(ax,X,Y,MarkerSize,markerColor,'filled','MarkerFaceAlpha',markerAlpha)
    scatter(ax,X,Y,MarkerSize,markerColor,'filled','MarkerFaceAlpha',markerAlpha); hold on;
    colormap(ax,"jet")
    %title(ax,varNameForSpots,'FontSize',8,'Interpreter','none')
    if uniformMaxX
        ax.XLim=[-maxX/50 maxX];
    else
        ax.XLim(1)=-ax.XLim(2)/50;
    end

    if uniformMaxY
        ax.YLim=[-maxY/50 maxY];
    else
        %ax.YLim=[-prctile(Y,yPercentile)/50 prctile(Y,yPercentile)*maxYPercentileFactor];
        ax.YLim(1)=-ax.YLim(2)/50;
    end
    if mean(diff(ax.YAxis.TickValues))<0.1
        ax.YAxis.TickLabelFormat='%.2f';
    else
        ax.YAxis.TickLabelFormat='%.1f';
    end
    if any(X>ax.XLim(2))
        fprintf('WARNING for %14s: %i of %i points are >maximum X shown (%i) with mean=%.1f\n',varNameForSpots,sum(any(X>ax.XLim(2))),length(X),ax.XLim(2),mean(X(X>ax.XLim(2))))
    end
    if any(Y>ax.YLim(2))
        fprintf('WARNING for %14s: %i of %i points are >maximum Y shown (%i) with mean=%.2f\n',varNameForSpots,sum(any(Y>ax.YLim(2))),length(Y),ax.YLim(2),mean(Y(Y>ax.YLim(2))))
    end

    % add text description of dataset to upper right or title
    rowTchannelInfo=find(strcmp(TchannelInfo.channelName,varNameForSpots));
    if textInAxis
        xText=ax.XLim(2)*0.99;
        yText=ax.YLim(2)*0.98;
        ht=text(xText,yText,           [TchannelInfo.gene{rowTchannelInfo},', ',TchannelInfo.amplifierSet{rowTchannelInfo}]);
        ht.HorizontalAlignment='right'; ht.FontSize=8;
        ht=text(xText,yText-textYspacing*diff(ax.YLim),      [TchannelInfo.readoutDye{rowTchannelInfo},', ',TchannelInfo.exposureTime{rowTchannelInfo}]);
        ht.HorizontalAlignment='right'; ht.FontSize=8;
        % show R^2
        RsquaredText=sprintf('R^{2}=%.3f',LM.Rsquared.Ordinary);
        ht=text(xText,yText-2*textYspacing*diff(ax.YLim),RsquaredText);
        ht.HorizontalAlignment='right'; ht.FontSize=8;

    else % put it in title
        %titleStr=sprintf('%s, %s, %s, %s\nR^{2}=%.3f',TchannelInfo.gene{rowTchannelInfo},TchannelInfo.amplifierSet{rowTchannelInfo},TchannelInfo.readoutDye{rowTchannelInfo},TchannelInfo.exposureTime{rowTchannelInfo},LM.Rsquared.Ordinary);
        title(['{\it',TchannelInfo.gene{rowTchannelInfo},'}, ',TchannelInfo.amplifierSet{rowTchannelInfo},', ',TchannelInfo.readoutDye{rowTchannelInfo},', ',TchannelInfo.exposureTime{rowTchannelInfo}],'FontSize',8,'FontWeight','normal')
        %text(0.3*ax.XLim(2),ax.YLim(2)*.98,sprintf('\nR^{2}=%.3f',LM.Rsquared.Ordinary),'FontSize',8)
        ht=text(0.98*ax.XLim(2),0,sprintf('\nR^{2}=%.3f',LM.Rsquared.Ordinary),'FontSize',8);
        ht.HorizontalAlignment='right'; ht.VerticalAlignment='bottom'; ht.FontSize=8;        
    end
    if plotWithRegressionOnAllData
        lmCoeff=LM.Coefficients.Estimate; % [intercenpt;slope]
        linRegress=@(x) lmCoeff(1) + lmCoeff(2)*x;
        fplot(linRegress,'-','Color','g')
    end

    % colorbar
    if strcmp(varNameForSpots,'R3_CY3_250ms_ITGA3')
        hColorbar=colorbar; % same colorbar for all of them (0 to 1)
        %hColorbar=findobj(t.Children,'Type','Colorbar');
        hColorbar.FontSize=8;
        hColorbar.Location="north";
        hColorbar.Position(1:4)=[0.4 0.99 0.2 .01];
        hColorbar.Ticks=[0:0.2:1];
        hColorbar.Label.String=densityColorLabel;
        hColorbar.Label.FontSize=8;
        hColorbar.Label.Position(2)=-2;
    end
end

if uniformMaxX
   plotName=sprintf('SpotVsMeanIntensity1_MaxX=%i',maxX);
else
    plotName='SpotVsMeanIntensity1';
end
if uniformMaxY
    plotName=sprintf('%s_MaxY=%0.2f',plotName,maxY);
end

writetable(TdataoutPlot,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)

exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']),'ContentType','vector')
exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.tif']),'ContentType','image','Resolution',300)



%% Plot 2: Spot count vs. Signal - PrevStripSignal
tic
TchannelsForThisPlot=TchannelsForPlot2;
f=figure(2);clf;
f.Position=FigurePosition;
numPlots=height(TchannelsForThisPlot);
numPlotCols=3;
numPlotRows=4;
t=tiledlayout(numPlotRows,numPlotCols);
%t=tiledlayout(1,1);
t.TileIndexing='columnmajor';
t.TileSpacing='compact';
t.Padding='normal';
xlabel(t,"mean background-subtracted fluorescence minus previous post-strip cycle's background-subtracted fluorescence (AU)",'FontSize',8)
ylabel(t,'spot count per \mum^{2}','FontSize',8)

% blank first few columns
for iPlotBlank=1:4 % fist column is blank
     ax=nexttile(t);
     ax.Visible='off';
end

% for simple plot data output
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig12_CountVsSig_wSubtract.xlsx';
outdataSheet='SupFig12_CountVsSig_wSubtract';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.
TdataoutPlot=Tcell(:,1:9);

fprintf('------------ plot 2 ------------\n')
for iPlot=1:numPlots
    %ax=nexttile(t,iPlot); % this weirdly has a bug when specifying the tile with columnmajor indexing, see https://www.mathworks.com/matlabcentral/discussions/highlights/134368-new-in-r2021a-improvements-to-tiled-chart-layout/40538
    ax=nexttile(t);
    ax.FontSize=8; ax.XLabel.FontSize=8; ax.YLabel.FontSize=8; hold on;
        
    fprintf('plotting iPlot=%i\n',iPlot)
    varNameForSpots=TchannelsForThisPlot.channelName{iPlot};
    prevStripChannel=TchannelsForThisPlot.prevStripChannel{iPlot};
    varNameForMeanThisCycle=['sig_',varNameForSpots];
    varNameForMeanPrevStripCycle=['sig_',prevStripChannel];


    X=Tcell.(varNameForMeanThisCycle) - Tcell.(varNameForMeanPrevStripCycle);
    Y=Tcell.(varNameForSpots)./Tcell.Area_Um2;

    if plotSparse
        % downsample the points by sparsityFactor
        rng(0);
        indSparse=sort(randperm(length(X),round(length(X)/sparsityFactor)))';
        X=X(indSparse);
        Y=Y(indSparse);
    end
    LM=fitlm(X,Y);

    TdataoutPlot.([varNameForSpots,'_mean'])=X;
    TdataoutPlot.([varNameForSpots,'_spot'])=Y;

    if plotWithDensity
        %Vq=ksdensity([X,Y],[X,Y]); % SLOW - use pts to evaluate kernel density are the same locations as the points themselves
        % faster way to use ksdensity is to use meshgrid of points and later interpolate density at the x,y points
        gridx=unique(prctile(X,0:1:100));
        gridy=unique(prctile(Y,0:1:100));
        [Xi,Yi]=meshgrid(gridx,gridy);
        [fdens,~]=ksdensity([X,Y],[Xi(:),Yi(:)]);
        Vq=interp2(Xi,Yi,reshape(fdens,size(Xi)),X,Y);%interpolate at the points of X,Y
        markerColor=Vq/max(Vq); % 1 is highest density
    else
        markerColor=markerColorUniform;
    end
    %scatter(ax,X,Y,MarkerSize,markerColor,'filled','MarkerFaceAlpha',markerAlpha)
    scatter(ax,X,Y,MarkerSize,markerColor,'filled','MarkerFaceAlpha',markerAlpha); hold on;
    colormap(ax,"jet")
    %title(ax,varNameForSpots,'FontSize',8,'Interpreter','none')
    if uniformMaxX
        ax.XLim=[-maxX/50 maxX];
    else
        ax.XLim(1)=-ax.XLim(2)/50;
    end

    if uniformMaxY
        ax.YLim=[-maxY/50 maxY];
    else
        %ax.YLim=[-prctile(Y,yPercentile)/50 prctile(Y,yPercentile)*maxYPercentileFactor];
        ax.YLim(1)=-ax.YLim(2)/50;
    end
    if mean(diff(ax.YAxis.TickValues))<0.1
        ax.YAxis.TickLabelFormat='%.2f';
    else
        ax.YAxis.TickLabelFormat='%.1f';
    end

    if any(X>ax.XLim(2))
        fprintf('WARNING for %14s: %i of %i points are >maximum X shown (%i) with mean=%.1f\n',varNameForSpots,sum(any(X>ax.XLim(2))),length(X),ax.XLim(2),mean(X(X>ax.XLim(2))))
    end
    if any(Y>ax.YLim(2))
        fprintf('WARNING for %14s: %i of %i points are >maximum Y shown (%i) with mean=%.2f\n',varNameForSpots,sum(any(Y>ax.YLim(2))),length(Y),ax.YLim(2),mean(Y(Y>ax.YLim(2))))
    end

     % add text description of dataset to upper right or title
    rowTchannelInfo=find(strcmp(TchannelInfo.channelName,varNameForSpots));
    if textInAxis
        xText=ax.XLim(2)*0.99;
        yText=ax.YLim(2)*0.98;
        ht=text(xText,yText,           [TchannelInfo.gene{rowTchannelInfo},', ',TchannelInfo.amplifierSet{rowTchannelInfo}]);
        ht.HorizontalAlignment='right'; ht.FontSize=8;
        ht=text(xText,yText-textYspacing*diff(ax.YLim),      [TchannelInfo.readoutDye{rowTchannelInfo},', ',TchannelInfo.exposureTime{rowTchannelInfo}]);
        ht.HorizontalAlignment='right'; ht.FontSize=8;
        % show R^2
        RsquaredText=sprintf('R^{2}=%.3f',LM.Rsquared.Ordinary);
        ht=text(xText,yText-2*textYspacing*diff(ax.YLim),RsquaredText);
        ht.HorizontalAlignment='right'; ht.FontSize=8;

    else % put it in title
        %titleStr=sprintf('%s, %s, %s, %s\nR^{2}=%.3f',TchannelInfo.gene{rowTchannelInfo},TchannelInfo.amplifierSet{rowTchannelInfo},TchannelInfo.readoutDye{rowTchannelInfo},TchannelInfo.exposureTime{rowTchannelInfo},LM.Rsquared.Ordinary);
        title(['{\it',TchannelInfo.gene{rowTchannelInfo},'}, ',TchannelInfo.amplifierSet{rowTchannelInfo},', ',TchannelInfo.readoutDye{rowTchannelInfo},', ',TchannelInfo.exposureTime{rowTchannelInfo}],'FontSize',8,'FontWeight','normal')
        %text(0.4*ax.XLim(2),ax.YLim(2)*.98,sprintf('\nR^{2}=%.3f',LM.Rsquared.Ordinary),'FontSize',8)
        ht=text(0.98*ax.XLim(2),0,sprintf('\nR^{2}=%.3f',LM.Rsquared.Ordinary),'FontSize',8);
        ht.HorizontalAlignment='right'; ht.VerticalAlignment='bottom'; ht.FontSize=8;        
    end


    if plotWithRegressionOnAllData
        lmCoeff=LM.Coefficients.Estimate; % [intercenpt;slope]
        linRegress=@(x) lmCoeff(1) + lmCoeff(2)*x;
        fplot(linRegress,'-','Color','g')
    end

    % colorbar
    %if strcmp(varNameForSpots,'R3_CY3_250ms_ITGA3')
    if iPlot==1
        hColorbar=colorbar; % same colorbar for all of them (0 to 1)
        %hColorbar=findobj(t.Children,'Type','Colorbar');
        hColorbar.FontSize=8;
        hColorbar.Location="north";
        hColorbar.Position(1:4)=[0.4 0.99 0.2 .01];
        hColorbar.Ticks=[0:0.2:1];
        hColorbar.Label.String=densityColorLabel;
        hColorbar.Label.FontSize=8;
        hColorbar.Label.Position(2)=-2;
    end
end
toc

if uniformMaxX
   plotName=sprintf('SpotVsMeanIntensity2_PrevStripSubtracted_MaxX=%i',maxX);
else
    plotName='SpotVsMeanIntensity2_PrevStripSubtracted';
end
if uniformMaxY
    plotName=sprintf('%s_MaxY=%0.2f',plotName,maxY);
end

writetable(TdataoutPlot,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)

exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']),'ContentType','vector')
exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.tif']),'ContentType','image','Resolution',300)
