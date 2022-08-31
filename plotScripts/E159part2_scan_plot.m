% E159part2_scan_plot (QC plots)
% high-throughput scan

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
Tcond=readtable(fullfile(parentDir,'experiments/E159part2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'experiments/E159part2_scan_Tscan.xlsx'));
extractedDataDir=fullfile(parentDir,'extractedData/E159part2_scan',filesep);
plotDir=fullfile(parentDir,filesep,'plots/scan',filesep);

cd(parentDir)


%% plot QC choices: Excluding certain subregions
TcellRaw=readtable(fullfile(extractedDataDir,'TcellAll1.csv'));
minMeanUBC=25;

numScan=1;
fontSize=16;

f=figure(1);f.Position=[500 500 1500 500];
f.Name='QC mean UBC per subregion';
t=tiledlayout(f,numScan,3);  t.TileSpacing = 'tight'; t.Padding = 'normal';
title(t,sprintf('mean UBC per cell, with outline around mean >=%.2g for all rounds',minMeanUBC),'FontSize',fontSize,'FontWeight','bold')
xlabel(t,'Subregion Array Column','FontSize',fontSize)
ylabel(t,'Subregion Array Row','FontSize',fontSize)

iAx=0;
for scanID=1:numScan
    
    %Tsr=groupsummary(TcellRaw(TcellRaw.scanID==scanID,{'subregionInd','subregionRow','subregionCol','scanID','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC','passMinUBC','passPctUBC','passAllQC'}),'uniqueSubregionInd','mean');
    Tsr=groupsummary(TcellRaw(TcellRaw.scanID==scanID,{'subregionInd','subregionRow','subregionCol','scanID','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC'}),'subregionInd','mean');
    Tsr=renamevars(Tsr,{'mean_subregionRow','mean_subregionCol','mean_scanID'},{'subregionRow','subregionCol','scanID'});
    Tsr=join(Tsr,unique(TcellRaw(:,{'scanID','scanWell'}),'rows'));
    Tsr=movevars(Tsr,'scanWell','After','scanID');
    assert(length(unique(Tsr.scanWell))==1)
    scanWell=char(unique(Tsr.scanWell));
    
    matSz=[max(Tsr.subregionRow),max(Tsr.subregionCol)];
    
    for thisRound=1:3
        Rmat=nan(matSz(1),matSz(2));
        varStr=['mean_R',num2str(thisRound),'_YFP_UBC'];
        ind=sub2ind(matSz,Tsr.subregionRow,Tsr.subregionCol);
        Rmat(ind)=Tsr.(varStr);
        RmatAll(:,:,thisRound)=Rmat;
    end
    
    % passQC
    
    passQCImg=all(RmatAll>minMeanUBC,3);
    [B,L] = bwboundaries(passQCImg);
    assert(length(B)==1)
    XY=fliplr(B{1});
    XY2=reshape(repmat(XY',4,1),2,[])' + repmat([-0.5 -0.5; -0.5 0.5; 0.5 -0.5;0.5 0.5],size(XY,1),1);
    XY2=round(2*XY2)/2;
    k=boundary(XY2,1);
    XY3=XY2(k,:);
    polyKeep=polyshape(XY3);
    
    
    % colormap
    cmap=turbo(ceil(max(RmatAll,[],'all')));
    % plot
    for  thisRound=1:3
        iAx=iAx+1;
        ax(iAx)=nexttile(iAx); axis equal; axis square; hold on; set(ax(iAx),'CLim',[0 size(cmap,1)])
        axes(ax(iAx));h1=imshow(RmatAll(:,:,thisRound),[0 size(cmap,1)],'Colormap',cmap);
        set(ax(iAx),'XLim',[0.5 matSz(1)+0.5])
        set(ax(iAx),'YLim',[0.5 matSz(2)+0.5])
        set(ax(iAx),'YDir','reverse')
        p1=plot(polyKeep,'FaceAlpha',0,'EdgeAlpha',1,'LineWidth',2);
        %set(ax(iAx),'PositionConstraint','outerposition')
        title(sprintf('well %s, Round %i',scanWell,thisRound),'FontSize',fontSize,'FontWeight','normal')
    end
    
    cb=colorbar('FontSize',fontSize);
    cb.Label.String="mean UBC per cell"; %cb.Label.Rotation=-90;
end

if ~isdir(plotDir)
    mkdir(plotDir)
end
exportgraphics(f, fullfile(plotDir,[f.Name,'.eps']), 'ContentType', 'vector');

%% histogram of cell area

%% Plots for QC decisions: UBC counts vs. cellarea 
TcellRaw=readtable(fullfile(extractedDataDir,'TcellAll1_withQC.csv'));

minUBCPerUm2=0.025;

zoomIn=true; % turn to false for full (uncropped)
xMaxShowZoom=2000;
yMaxShowZoom=400;

fontSize=16;

if zoomIn==false
    figNum=3;
    markerAlpha=1;
else
    figNum=4;
    markerAlpha=0.5;
end

markerSize=5;

if zoomIn==false
    fName='UBC counts vs cell area';
else
    fName='UBC counts vs cell area (zoomed view)';
end

f = figure(figNum); f.Name=fName; clf; f.Position=[500 500 1500 500]; set(f,'visible','off'); ax(iAx)=gca;
t=tiledlayout(f,numScan,3);  t.TileSpacing = 'tight'; t.Padding = 'normal';

if zoomIn==false
    titleStr=sprintf('UBC counts (showing cutoff %.3f UBC per Um^2)',minUBCPerUm2);
else
    titleStr=sprintf('UBC counts (showing cutoff %.3f UBC per Um^2), zoomed view: Area [0 %.0f], UBC count [0 %.0f]',minUBCPerUm2,xMaxShowZoom,yMaxShowZoom);
end

title(t,titleStr,'FontSize',fontSize,'FontWeight','bold')

xlabel(t,'Cell area, Um^2','FontSize',fontSize)
ylabel(t,'UBC count','FontSize',fontSize)


idxLogical=startsWith(TcellRaw.Properties.VariableNames,'pass');
TcellRaw=convertvars(TcellRaw,TcellRaw.Properties.VariableNames(idxLogical),'logical');

Tcell=TcellRaw(TcellRaw.passSubregion,:);
%Tcell=TcellRaw;



xVar='Area_Um2';
yVarArray={'R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC'};
for iAx=1:3
yVar=yVarArray{iAx};

X=Tcell.(xVar);
Y=Tcell.(yVar);
C=Tcell.subregionCol;

% colormap
cmap=turbo(ceil(max(RmatAll,[],'all')));
ax(iAx)=nexttile(iAx); set(ax(iAx),'CLim',[0 size(cmap,1)])
sh=scatter(X,Y,markerSize,C,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',markerAlpha);
title(sprintf('Round %i',iAx),'FontWeight','normal');
hold on; %axis equal; axis square;

end

set(ax,'FontSize',fontSize)

cb=colorbar('FontSize',fontSize,'Limits',[1 max(C)]);
cb.Label.String="Subregion Array Column"; %cb.Label.Rotation=-90;




if zoomIn
    set(ax,'XLim',[0 xMaxShowZoom])
    set(ax,'YLim',[0 yMaxShowZoom])
else
    yLimMax=max(cellfun(@(x) x(2), get(ax,'YLim')));
    set(ax,'YLim',[0 yLimMax])
end

% plot minimum line
xLine=get(ax(1),'XLim');
yLine=xLine*minUBCPerUm2;
for iAx=1:3
plot(ax(iAx),xLine,yLine,'k')
end
set(ax,'Colormap',jet(max(C)))
linkaxes(ax)
set(f,'visible','on')

exportgraphics(f, fullfile(plotDir,[f.Name,'.eps']), 'ContentType', 'vector');

%% Plots for QC decisions: UBC/cellarea 
fontSize=16;

TcellRaw=readtable(fullfile(extractedDataDir,'TcellAll1_withQC.csv'));
idxLogical=startsWith(TcellRaw.Properties.VariableNames,'pass');
TcellRaw=convertvars(TcellRaw,TcellRaw.Properties.VariableNames(idxLogical),'logical');
Tcell=TcellRaw(TcellRaw.passSubregion,:);
markerAlpha=1;
markerSize=5;

xVar='R1_YFP_UBC';
yVar='R3_YFP_UBC';
X=Tcell.(xVar)./Tcell.Area_Um2;
Y=Tcell.(yVar)./Tcell.Area_Um2;
C=Tcell.subregionCol;
f = figure(2); f.Name='UBC vs per cell area.eps'; clf; set(f,'visible','off'); ax=gca;
scatter(X,Y,markerSize,C,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',markerAlpha);
hold on; axis equal; axis square;
cb=colorbar('FontSize',fontSize);
cb.Label.String="Subregion Array Column"; %cb.Label.Rotation=-90;
maxVal=max([[ax.XLim],[ax.YLim]]);
set(ax,'XLim',[0 maxVal])
set(ax,'YLim',[0 maxVal])

set(ax,'FontSize',fontSize)
plot(ax,ax.XLim,ax.XLim,'k')
text(ax,ax.XLim(2)*0.75,ax.XLim(2)*0.70,'y=x','FontSize',fontSize)
set(ax,'Colormap',jet(10))
set(ax,'YTick',get(ax,'XTick'))
xlabel(sprintf('%s counts per Um^2',xVar),'Interpreter','none');ylabel(sprintf('%s counts per Um^2',yVar),'Interpreter','none')
set(f,'visible','on')

exportgraphics(f, fullfile(plotDir,f.Name), 'ContentType', 'vector');
%% plot UBC correlations
afterQC=false;
zoomIn=true;
maxValShow=300;

fontSize=16;

if afterQC
    titleStr='UBC correlation between rounds (after QC filters)';
    Tcell=readtable(fullfile(extractedDataDir,'TcellAll2.csv'));
    figNum=6;
else % before QC
    titleStr='UBC correlation between rounds (before QC filters)';
    figNum=5;
    Tcell=readtable(fullfile(extractedDataDir,'TcellAll1.csv'));
end

if zoomIn==true
    titleStr=[titleStr,sprintf(', Zoom to UBC count [0 %.0f]',maxValShow)];
end

fName=titleStr;
f = figure(figNum); f.Name=fName; f.Position=[500 500 1500 500]; set(f,'visible','off');

C=Tcell.subregionCol; % color index

markerAlpha=0.5;
markerSize=5;
t=tiledlayout(1,3);  t.TileSpacing = 'tight';t.Padding = 'tight';
title(t,titleStr,'FontSize',fontSize,'FontWeight','bold')

ax(1)=nexttile(1); axis equal; axis square; hold on;
scatter(Tcell.R1_YFP_UBC,Tcell.R2_YFP_UBC,markerSize,C,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',markerAlpha)
xlabel('R1'); ylabel('R2')
title('UBC: R2 vs. R1','FontSize',fontSize)

ax(2)=nexttile(2); axis equal; axis square; hold on;
scatter(Tcell.R1_YFP_UBC,Tcell.R3_YFP_UBC,markerSize,C,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',markerAlpha)
xlabel('R1'); ylabel('R3')
title('UBC: R3 vs. R1','FontSize',fontSize)

ax(3)=nexttile(3); axis equal; axis square; hold on;
scatter(Tcell.R2_YFP_UBC,Tcell.R3_YFP_UBC,markerSize,C,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',markerAlpha)
xlabel('R2'); ylabel('R3')
title('UBC: R3 vs. R2','FontSize',fontSize)

% set xlim and ylim
maxVal=max([[ax.XLim],[ax.YLim]]);
if zoomIn==true
    set(ax,'XLim',[0 maxValShow])
    set(ax,'YLim',[0 maxValShow])
else
    set(ax,'XLim',[0 maxVal])
    set(ax,'YLim',[0 maxVal])
end

% draw X=Y
for i=1:length(ax)
    plot(ax(i),ax(i).XLim,ax(i).XLim,'k')
    text(ax(i),ax(i).XLim(2)*0.75,ax(i).XLim(2)*0.70,'y=x','FontSize',fontSize)
    set(ax(i),'Colormap',jet(10))
    %set(ax(i),'Colormap',jet(maxUniqueSubregionInd))
    set(ax(i),'FontSize',fontSize)
end

cb=colorbar('FontSize',fontSize);
cb.Label.String="Subregion Array Column"; %cb.Label.Rotation=-90;

set(f,'visible','on')

exportgraphics(f, fullfile(plotDir,[f.Name,'.eps']), 'ContentType', 'vector');

%% plot histograms of all RNA counts
fontSize=16;
plotLog=false;
FaceColor=[0, 150, 255]/255;

plotCountThresh=true;


Tcell=readtable(fullfile(extractedDataDir,'TcellAll2.csv'));
fishChannels=   {'R1_YFP_UBC','R1_CY3_WNT5A','R1_A594_DDX58','R1_CY5_AXL','R2_YFP_UBC','R2_CY3_NGFR','R2_A594_FN1','R2_CY5_EGFR','R3_YFP_UBC','R3_CY3_ITGA3','R3_A594_MMP1','R3_CY5_MITF'};
countThresholds=[nan           15             10             25           nan           30           100            5           nan           50            40             nan];
temp=split(fishChannels','_');
temp=temp(:,2)';
chanColorInd=cellfun(@(x) find(strcmp(x,{'YFP','CY3','A594','CY5'})) , temp);
assert(isequal(chanColorInd,repmat(1:4,1,3)));




titleStr=sprintf('Spot count histograms for %i cells',height(Tcell));
if plotLog
    figNum=7;
    titleStr=[titleStr,' (log plot)'];
else
    figNum=8;
    titleStr=[titleStr,' (linear plot)'];
end
fName=titleStr;

f = figure(figNum); f.Name=fName; f.Position=[1000 200 900 700]; set(f,'visible','on');
t=tiledlayout(4,3);  t.TileSpacing = 'tight';t.Padding = 'tight';
title(t,titleStr,'FontSize',fontSize,'FontWeight','bold')

iAx=0;
for row=1:4
    %channelIndThisRow=[(3*(row-1)+1):(3*(row))];
    channelIndThisRow=find(chanColorInd==row);
    for col=1:3
        iAx=iAx+1;
        ax(iAx)=nexttile(iAx); hold on;
        set(ax(iAx),'FontSize',fontSize)
        
        channelInd=channelIndThisRow(col);
        channelName=fishChannels{channelInd};
        
        RNAcounts=Tcell.(channelName);
        h1=histogram(ax(iAx),RNAcounts,60,'FaceColor',FaceColor,'EdgeColor','none','FaceAlpha',1);
        set(ax(iAx),'PositionConstraint','innerposition')
        
        if plotLog
           set(ax(iAx),'YScale','log') 
        end
        
        xlm=get(ax(iAx),'XLim');
        ylm=get(ax(iAx),'YLim');
        
        %th=title(channelName,'FontSize',fontSize,'FontWeight','normal','Interpreter','none','Position',titlePosition);
                    if plotLog
                        yTitlePos=10^(0.85*log10(ylm(2)));
                    else
                        yTitlePos=0.85*ylm(2);
                    end
                    titlePosition=[xlm(1)+ diff(xlm)*0.5, yTitlePos];
        tith=title(channelName,'FontSize',fontSize,'FontWeight','normal','Interpreter','none','Position',titlePosition);
        
        % plot threshold line
        if plotCountThresh
                xlm=get(ax(iAx),'XLim');
                ylm=get(ax(iAx),'YLim');
                thresh=countThresholds(channelInd);
                if ~isnan(thresh)
                    plot(ax(iAx),repmat(thresh,1,2),ylm,'k','LineWidth',2)
                    if plotLog
                        yLabelPos=10^(0.65*log10(ylm(2)));
                    else
                        yLabelPos=0.65*ylm(2);
                    end
                    xLabelPos=thresh+0.01*diff(xlm);
                    pctAbove=100*sum(RNAcounts>=thresh)/size(RNAcounts,1);
                    th=text(xLabelPos,yLabelPos,sprintf('%.0f  (%.2f%% above)',thresh,pctAbove),'FontSize',fontSize);
                end
        end
    end % end col loop
    
    if row==1
        % set UBC limits equal to eachother
        iAxToSet=(iAx-2):iAx;
        xLimMin=min(cellfun(@(x) x(1), get(ax(iAxToSet),'XLim')));
        xLimMax=max(cellfun(@(x) x(2), get(ax(iAxToSet),'XLim')));
        yLimMin=min(cellfun(@(x) x(1), get(ax(iAxToSet),'YLim')));
        yLimMax=max(cellfun(@(x) x(2), get(ax(iAxToSet),'YLim')));
        set(ax(iAxToSet),'XLim',[xLimMin xLimMax])
        set(ax(iAxToSet),'YLim',[yLimMin yLimMax])
    end
    
end
set(f,'visible','on');

exportgraphics(f, fullfile(plotDir,[f.Name,'.eps']), 'ContentType', 'vector');



















