% E157part2BothRep_VsSmFISH_plot
% clampFISH 2.0 vs. smFISH correlation for rep1 and rep2
%
% note: E157part2rep2, 10X data at 1x1 binning was unfortunately not taken before
% readout stripping and smFISH for DDX58. Only have after stripping and smFISH and
% subsequent readout, which is different than rest.
%
% Change Ymag as follows:
%   Ymag='20X'; %--> For main figure
%   Ymag='10X'; %--> For supplement. 

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)

Xmag='60X';
Xmethod='conventional single-molecule RNA FISH';
Ymag='20X'; % Ymag='20X'--> For main figure; Ymag='10X' --> For Extended data
Ymethod='clampFISH 2.0';


% make Tplots, based on Tcond
% Rep1: Tcond --> Tplots
extractedDataTable='paper/extractedData/E157_validation/VsSmFISH/TcellMapBetweenMagnifications2.csv';
Tcond_rep1=readtable('paper/experiments/E157_validation_Tcond.xlsx');
Tcond_rep1.rep=repmat(1,height(Tcond_rep1),1);
Tcond_rep1.extractedDataTable=repmat({extractedDataTable},height(Tcond_rep1),1);
rep1_condIDs_toPlot=[12,14,15];
plotRows=           [1 1 1];
plotCols=           [1 2 3];
assert(isequal(Tcond_rep1.condID,(1:height(Tcond_rep1))')) % rows same as condID
Tplots_rep1=[array2table([plotRows',plotCols'],'VariableNames',{'plotRow','plotCol'}),...
    Tcond_rep1(rep1_condIDs_toPlot,{'rep','condID','cellType','AmpRound','ReadoutSetName','extractedDataTable','spots_10X_channel','spots_10X_threshold','spots_20X_channel','spots_20X_threshold','spots_60X_channel','spots_60X_threshold'})];

% Rep2: Tcond --> Tplots
extractedDataTable='paper/extractedData/E157part2rep2_VsSmFISH/TcellMapBetweenMagnifications2.csv';
Tcond_rep2=readtable('paper/experiments/E157part2rep2_VsSmFISH_Tcond.xlsx');
Tcond_rep2.rep=repmat(2,height(Tcond_rep2),1);
Tcond_rep2.extractedDataTable=repmat({extractedDataTable},height(Tcond_rep2),1);
rep2_condIDs_toPlot=[3,2,1];
plotRows=           [2 2 2];
plotCols=           [1 2 3];
assert(isequal(Tcond_rep2.condID,(1:height(Tcond_rep2))')) % rows same as condID
Tplots_rep2=[array2table([plotRows',plotCols'],'VariableNames',{'plotRow','plotCol'}),...
    Tcond_rep2(rep2_condIDs_toPlot,{'rep','condID','cellType','AmpRound','ReadoutSetName','extractedDataTable','spots_10X_channel','spots_10X_threshold','spots_20X_channel','spots_20X_threshold','spots_60X_channel','spots_60X_threshold'})];
% combine into Tplots
Tplots=[Tplots_rep1;Tplots_rep2];
Tplots.plotID=(1:height(Tplots))';
Tplots=movevars(Tplots,'plotID','Before',1);
Tplots

% Export simple data used for plot
outdataFile='SourceData_Fig2.xlsx';
outdataSheet='SourceData_Fig2';
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);

% make TcellData from TcellMapBetweenMagnifications2.csv files
TcellData=table();
for i=1:height(Tplots)
    plotID=Tplots.plotID(i);
    rep=Tplots.rep(i);
    condID=Tplots.condID(i);
    
    extractedDataTable=Tplots.extractedDataTable{i};
    TcellData_Temp=readtable(extractedDataTable);
    TcellData_Temp=TcellData_Temp(TcellData_Temp.condID==condID,:);
    
    TcellData_Temp.plotID=repmat(plotID,height(TcellData_Temp),1);
    TcellData_Temp.rep=repmat(rep,height(TcellData_Temp),1);
    TcellData_Temp.condID=repmat(condID,height(TcellData_Temp),1);
    TcellData_Temp=movevars(TcellData_Temp,{'plotID','rep','condID'},'Before',1);
    TcellData=[TcellData;TcellData_Temp];
end
TcellData=join(TcellData,Tplots(:,{'rep','condID','ReadoutSetName'}),'Keys',{'rep','condID'});
TcellData=movevars(TcellData,'ReadoutSetName','After','condID');
%TcellData

% plotting inputs
plotOutDir='paper/plots/E157part2BothRep_VsSmFISH';


plotOutName=['E157_VsSmFISH',Ymag];

% Export simple data used for plot
writetable(TcellData,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)


channelLastLetterMap=containers.Map({'','a','b','c','d','e'},{25,50,100,250,500,1000});

%XorYMaxMode='auto';
XorYMaxMode='manual';
maxXorY_manual=[1650 949 30 1650 949 30];

% check max data
condIDsToCheck=Tplots.condID(strcmp(Tplots.ReadoutSetName,'EGFR'));
idx=ismember(TcellData.condID,condIDsToCheck);
maxXorY_check=max([TcellData.numSpots_10X(idx),TcellData.numSpots_20X(idx),TcellData.numSpots_60X(idx)],[],'all')


figNum=12

addTextboxWithExposureTime=true;

jitterMax=0.2;

fontSize=8;
regressionFontSize=8;
markerSize=3;
markerColor=[0 0 250]/255;

markerAlpha=0.2;

regressionTextInTopLeft=true;

plot_regression=true;

fh=figure(figNum); clf;
set(fh,'Position',[1400 200 375 290],'Units','points')
t=tiledlayout(max(Tplots.plotRow),max(Tplots.plotCol));
t.TileSpacing='compact';
t.Padding='compact';

xlabel(t,sprintf('%s %s spot count',Xmag,Xmethod),'FontSize',fontSize,'FontWeight','normal')
ylabel(t,sprintf('%s %s spot count',Ymag,Ymethod),'FontSize',fontSize,'FontWeight','normal')

for i=1:height(Tplots)
    plotID=Tplots.plotID(i);
    plotRow=Tplots.plotRow(i);
    plotCol=Tplots.plotCol(i);
    
    rep=Tplots.rep(i);
    condID=Tplots.condID(i);
    geneName=Tplots.ReadoutSetName{i};
    AmpRound=Tplots.AmpRound{i};
    cellType=Tplots.cellType{i};
    
    ax=nexttile(t,plotID); hold on; axis equal;
    ax.FontSize=fontSize;
    
    % TITLE
    %if plotRow==1
    %hSgtitle=sgtitle(sprintf('%s %s, drug-%s cells',geneName,AmpRound,cellType,condID)); %verbose
    %set(hSgtitle,'FontSize',18)
    %set(hSgtitle,'FontWeight','bold')
    %end
    
    TcellData_ThisPlot=TcellData(TcellData.plotID==plotID,:);
    X=TcellData_ThisPlot.(['numSpots_',Xmag]);
    Y=TcellData_ThisPlot.(['numSpots_',Ymag]);
    
    % print out % sensitivity results
    fprintf('rep %i, %5s: clampFISH 2.0 detected %6i out of %6i smFISH spots (%.0f%%)\n',rep,geneName,sum(Y),sum(X),100*sum(Y)/sum(X))

    rng(0);
    Xjitter=(rand(size(X))-0.5)*2*jitterMax;
    Yjitter=(rand(size(X))-0.5)*2*jitterMax;

    scatter(X+Xjitter,Y+Yjitter,markerSize,markerColor,'filled','MarkerFaceAlpha',markerAlpha)%,'Size', markerSize)
    
    axis equal;
    
    switch XorYMaxMode
        case 'auto'
            maxXorY=max([ax(subplotCounter).XLim,ax(subplotCounter).YLim]);        
        case 'manual'
            maxXorY=maxXorY_manual(i);
    end
    
    ax.XLim=[0 maxXorY];
    ax.YLim=[0 maxXorY];
    ax.FontSize=fontSize;
    
    if diff(ax.XTick(1:2))>diff(ax.YTick(1:2))
        ax.YTick=ax.XTick;
    else
        ax.XTick=ax.YTick;
    end
    

    
    % TITLES
        channelNameY=Tplots.(['spots_',Ymag,'_channel']){i};
        thresholdY=Tplots.(['spots_',Ymag,'_threshold'])(i);
        
        YchannelLastLetter=replace(channelNameY,{'gfp','tmr','alexa','cy'},{''});
        YchannelExposureTime=channelLastLetterMap(YchannelLastLetter);
        
        if isnumeric(thresholdY)
            thresholdY=num2str(thresholdY);
        elseif iscell(thresholdY)
            thresholdY=thresholdY{:};
        end
        
        channelNameX=Tplots.(['spots_',Xmag,'_channel']){i};
        thresholdX=Tplots.(['spots_',Xmag,'_threshold'])(i);
        if isnumeric(thresholdX)
            thresholdX=num2str(thresholdX);
        elseif iscell(thresholdX)
            thresholdX=thresholdX{:};
        end
        
%             titleStr=sprintf('%s (%s, %ims, >=%s)',Ymag,channelNameY,YchannelExposureTime,thresholdY);
%             titleStrAddition=['{\it',sprintf('%s',geneName),'} ']; % gene name
%             titleStr=[titleStrAddition,titleStr];
            %hTitle=title(titleStr);
            %set(hTitle,'FontWeight','normal')
            if plotCol==2
                hTitle=title(sprintf('replicate %i',rep));
                set(hTitle,'FontWeight','normal')
            end
            
    
    % X labels
            %xlabel('conventional single-molecule RNA FISH spot count')
            %xlabel(sprintf('%s %s spot count',Xmag,Xmethod))    
    % Y label
%         if plotCol==1
%             ylabel(ax,sprintf('%s %s\nspot count',Ymag,Ymethod))
%         end
    
    % REGRESSION LINE
    
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
    
    %regressionText=sprintf('R^2=%.3f\nm= %.2f\nb= %.1f',r^2,m,b);
        if sign(b)==1 || b==0
            AddOrSubtractB='+';
        elseif sign(b)==-1
            AddOrSubtractB='-';
        else
            error('huh')
        end
        
        %regressionText=sprintf('R^2=%.3f\ny=%.2fx %s %.1f',Rsq,m,AddOrSubtractB,abs(b));
        regressionText=sprintf('y=%.2fx %s %.1f\nR^2=%.3f',m,AddOrSubtractB,abs(b),Rsq);
        geneStr=['{\it',sprintf('%s',geneName),'} ']; % gene name
        inPlotStr=sprintf('%s\n%s',geneStr,regressionText);
        text(ax.XLim(2)*.05,ax.YLim(2)*0.8,inPlotStr,'FontSize',regressionFontSize)

    
    
%     % add exposure time
%     if addTextboxWithExposureTime
%         text(ax,0.1*xlim(gca),0.6*ylim(gca),sprintf('%ims',YchannelExposureTime),'FontSize',fontSize)
%     end

        if condID==1 && strcmp(Ymag,'10X')
        % E157part2rep2, 10X data at 1x1 binning was unfortunately not taken before
        % readout stripping and smFISH for DDX58. Only have after stripping and smFISH and
        % subsequent readout, which is different than rest.
        ax.Visible='off';
        delete(ax.Children)
        end

end % end Tplots loop

% rep1
TcellData_thisRep=TcellData(TcellData.rep==1,:);
TcountStats=grpstats(TcellData(:,{'rep','ReadoutSetName','numSpots_60X','numSpots_20X'}),{'rep','ReadoutSetName'});
TcountStats=renamevars(TcountStats,'GroupCount','numCells')
rep1_aggregate_efficiencyPct=100*sum(TcountStats.mean_numSpots_20X(TcountStats.rep==1))/sum(TcountStats.mean_numSpots_60X(TcountStats.rep==1));
rep2_aggregate_efficiencyPct=100*sum(TcountStats.mean_numSpots_20X(TcountStats.rep==2))/sum(TcountStats.mean_numSpots_60X(TcountStats.rep==2));
fprintf('rep 1 aggregated efficiency: %.0f%%\n',rep1_aggregate_efficiencyPct)
fprintf('rep 2 aggregated efficiency: %.0f%%\n',rep2_aggregate_efficiencyPct)
% export
if ~isdir(plotOutDir)
    mkdir(plotOutDir)
end
exportgraphics(fh, fullfile(plotOutDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector'); % with spots

