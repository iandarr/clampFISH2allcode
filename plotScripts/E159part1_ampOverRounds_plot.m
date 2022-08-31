% E159part1_ampOverRounds_plot
% amplification to rounds 1,2,4,6,8,10
% also, amplification to rounds 2,8 without click reaction

% makes 3 plots:
%   1. For main fig: amplification over rounds and +/- click for UBC only
%   2. For supplementary fig: amplification over rounds for 4 genes
%   3. For supplementary fig: +/- click for 4 genes
%
% this script requires violinplot:
%           Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
%           https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847

%% Specify directories and import data
clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
parentDirForExport='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
cd(parentDir)
Tcond=readtable('paper/experiments/E159part1_ampclick_Tcond.xlsx');
% Import Data
cd(parentDir)
extractedDataDirForAllCond='paper/extractedData/E159_scan_amp/analysis1_amp';
plotDir=['paper/plots/E159part1_amplification_and_click'];

Tspots=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tspots4AllNormalized.csv'),'ReadVariableNames', true);
%Tcell=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tcell2All.csv'),'ReadVariableNames', true);

%% universal settings for all 3 plots
% regression line through medians (amplification plot only)
plot_medians_regression_line=true;
medians_textlabels=true;
fontSize=8;
regressionFontSize=8;

% font sizes
axisLabelFontSize=8;
titleFontSize=8;

% appearance of violins
ViolinAlpha=1;
clickPositiveColor=[36 255 36]/255;
clickNegativeColor=[0.8 0.8 0.8];

% appearanch of boxplots
removeWhiskers=true;

% appearance of median marker
MedianMarker='o'; %default 'o'
MedianLineWidth=0.5; %default 0.5
MedianSizeData=14; % default 36
MedianColor=[1 1 1]; %[1 1 1]=white

% appearance of interquartile range box
BoxWidth=0.02;
BoxColor=[0.4 0.4 0.4];


%% 1. Main figure plot: amplification and +/- click plots for UBC

numPlots=2;

% [amplification, click]
plotPositionNums=[1 2];

figureWidthPoints=180;
figureHeightPoints=180;

cd(parentDir)
figure(1);
fh=gcf; clf; hold on;
set(fh,'Position',[1500 200 figureWidthPoints figureHeightPoints],'Units','points')
delete(gca)

for iPlot=1:numPlots
    
    channelBaseName='gfp';
    positionNum=plotPositionNums(iPlot);
    
    if iPlot==1 %amplification
        plotCondIDs=[1 2 3 4 5 6];
    elseif iPlot==2 %click
        plotCondIDs=[2 7 5 8];
    else
        error('unknown case')
    end
    
    
    ampRoundsOrdered=Tcond.AmpRound(ismember(Tcond.condID,plotCondIDs));
    if iPlot==1
        assert(isequal(ampRoundsOrdered,[1 2 4 6 8 10]'));
    elseif iPlot==2
        assert(isequal(ampRoundsOrdered,[2 8 2 8]'));
    end
    
    indTspots=all([strcmp(Tspots.channelBaseName,channelBaseName),ismember(Tspots.condID,plotCondIDs)],2);
    if iPlot==1 %amplification
        X=Tspots.AmpRound(indTspots);
        assert(isequal(unique(X),[1 2 4 6 8 10]'));
        specialClicks={'','','','','',''}';
    elseif iPlot==2 %click
        X=Tspots.condName(indTspots);
        assert(isequal(unique(X),{'R2';'R2_NoCu';'R8';'R8_NoCu'}));
        specialClicks={'','NoCu','','NoCu'}';
    else
        error('unknown case')
    end
    Y=Tspots.intensitiesNorm(indTspots);
    
    % create subplots within the figure, 1x4 arrangement, and select one
    %subplot(1,2,positionNum)
    %ax(iPlot)=gca;
    %ax(iPlot).PositionConstraint = 'outerposition';
    
    
    % get position of new axes
    XleftSpace=0.11;
    XaxesWidthAmplification=0.565;
    XbetweenAxesSpace=0.05;
    XaxesWidthClick=XaxesWidthAmplification * 4.5/11;
    YbotSpace=0.15;
    YaxesHeight=0.87;
    if iPlot==1 %amplification
        XPosAxes= XleftSpace;
        AxesPosition=[XPosAxes YbotSpace XaxesWidthAmplification YaxesHeight];
    elseif iPlot==2 %click
        XPosAxes= XleftSpace + XaxesWidthAmplification + XbetweenAxesSpace;
        AxesPosition=[XPosAxes YbotSpace XaxesWidthClick YaxesHeight];
    end
    

    
    % create axes
    ax(iPlot)=axes('Position',AxesPosition,'PositionConstraint','innerposition');
    
    % XLim
    if iPlot==1
        ax(iPlot).XLim=[0 11];
    elseif iPlot==2
        ax(iPlot).XLim=[0 4.5];
    end
    
    %YLim
    ax(iPlot).YLim=[0.5 220]; % set equal to reference
    % log scale
    ax(iPlot).YScale='log';

    % plot violins
    vs=violinplot(Y,X,'ShowData',false,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    
    % recolor the violins
    for iViolin=1:length(vs)
        specialClick=specialClicks{iViolin};
        switch specialClick
            case ''
                vs(iViolin).ViolinPlot.FaceColor=clickPositiveColor;
            case 'NoCu'
                vs(iViolin).ViolinPlot.FaceColor=clickNegativeColor;
            otherwise
                error('case not handled')
        end
    end
    
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs)
            vs(iViolin).WhiskerPlot.LineStyle='none';
        end
    end
    
    % make violin X's be centered around....
    XDistanceToAddArray=nan(1,length(vs));
    if iPlot==1 %amplification
        % AmpRound (1,2,4,6,8,10)
        for iViolin=1:length(vs)
            thisAmpRound=ampRoundsOrdered(iViolin);
            XDistanceToAddArray(iViolin)=thisAmpRound - iViolin;
        end
        vs=XshiftViolinplot(vs,XDistanceToAddArray);
    elseif iPlot==2
        % center around 1.5 and 3.5. Currently centered around 1,2,3,4
        XDistanceToAddArray=[0.15 -0.15 0.15 -0.15];
        vs=XshiftViolinplot(vs,XDistanceToAddArray);
    end
    
    % median marker
    Xmedians_ViolinPlot=nan(length(plotCondIDs),1);
    Ymedians_ViolinPlot=nan(length(plotCondIDs),1);
    for iViolin=1:length(vs)
        vs(iViolin).MedianPlot.Marker=MedianMarker;
        vs(iViolin).MedianPlot.LineWidth=MedianLineWidth;
        vs(iViolin).MedianPlot.SizeData=MedianSizeData;
        Xmedians_ViolinPlot(iViolin)=vs(iViolin).MedianPlot.XData;
        Ymedians_ViolinPlot(iViolin)=vs(iViolin).MedianPlot.YData;
    end
    
    % medians for regression and labels
    Ymedians=nan(length(plotCondIDs),1);
    for iCondID=1:length(plotCondIDs)
        condID=plotCondIDs(iCondID);
        YdataThisGeneThisCond=Tspots.intensitiesNorm(all([strcmp(Tspots.channelBaseName,channelBaseName),Tspots.condID==condID],2));
        Ymedians(iCondID)=median(YdataThisGeneThisCond);
    end
    
    % check that medians in ViolinPlot are same as medians for regression
    assert(isequal(Ymedians_ViolinPlot,Ymedians))
    
    % plot medians regression line
    if (iPlot==1)&&plot_medians_regression_line
        XdataReg=ampRoundsOrdered(2:end);
        YdataReg=Ymedians(2:end);
        %f = fit(XdataReg,YdataReg,'power2');
        coeff=polyfit(XdataReg',log2(YdataReg'),1); % use data before the jitter
        m=coeff(1); b=coeff(2);
        %XregLim=[ampRoundsOrdered(2) ampRoundsOrdered(end)];
        XregLim=XdataReg;
        YregLim=2.^(m*XregLim+b);
        plot(ax(iPlot),XregLim,YregLim,'-','Color','k')
        regressionText=sprintf('y=%.2f\\times%.3f^{Round}',2^b,2^m);
        %text(ax(iPlot),.30,0.17,regressionText,'FontSize',regressionFontSize,'Interpreter','none','Units', 'normalized')
        text(ax(iPlot),.30,0.17,regressionText,'FontSize',regressionFontSize,'Units', 'normalized')
    end
    
    % medians text labels
    if medians_textlabels
        for iCondID=1:length(plotCondIDs)
            
            if (iPlot==2)&&(iCondID==1)
                thisMedianLabelOffsetX=-1.5;
            elseif (iPlot==2)&&(iCondID==3)
                thisMedianLabelOffsetX=-2;
            else
                thisMedianLabelOffsetX=0.35;
            end
            
            thisMedianLabelOffsetY= - Ymedians_ViolinPlot(iCondID)*0.05;
            
            medianLabelX=Xmedians_ViolinPlot(iCondID) + thisMedianLabelOffsetX;
            medianLabelY=Ymedians_ViolinPlot(iCondID) + thisMedianLabelOffsetY;
            

            text(ax(iPlot),medianLabelX,medianLabelY,num2str(Ymedians_ViolinPlot(iCondID),'%0.1f'),'FontSize',regressionFontSize)
        end
    end
    

    

    
    % X Tick Labels
    if iPlot==1
        ax(iPlot).XAxis.TickValues=ampRoundsOrdered;
    elseif iPlot==2
        ax(iPlot).XTick=[1.5 3.5];
        ax(iPlot).XTickLabel=[2 8];
    end
    
    % Y Tick Labels (turn off for right-side subplots)
    ax(iPlot).YTick=[1 2 4 8 16 32 64 128];
    if iPlot==2
        ax(iPlot).YTickLabel=[];
    end
    ax(iPlot).YTick=[1 2 4 8 16 32 64 128];
    ax(iPlot).YMinorTick='off';
    
    ax(iPlot).FontSize=fontSize;
   
    % make axis background invisible 
    set(ax(iPlot), 'color', 'none');
    
    % make y-axis invisible
    if iPlot==2
        ax(iPlot).YAxis.Visible = 'off';
    end
    
    % xlabel for whole figure
    if iPlot==1
        xlabelText='Round of amplification';
        text(0.4, -0.1, xlabelText, 'Units', 'normalized', 'FontSize',axisLabelFontSize,'FontWeight','normal')
    end
    
    % ylabel for whole figure
    if iPlot==1
        ylabelText= text(-0.15,0.2,'Spot intensity, normalized','Units','normalized','FontSize',axisLabelFontSize,'FontWeight','normal');
        set(ylabelText,'Rotation',90)
    end
    
    % legend
	if iPlot==2
        text(-0.08, 0.94,     '(+) Click','Color',clickPositiveColor,'FontSize',fontSize,'Units','normalized','FontWeight','bold')
        text(-0.05, 0.94-0.07,'(-) Click','Color',clickNegativeColor,'FontSize',fontSize,'Units','normalized','FontWeight','bold')
    end
end

cd(parentDirForExport)
if ~isfolder(plotDir)
    mkdir(plotDir)
end
exportgraphics(gcf, fullfile(plotDir,'E159_amplification_click_UBC.eps'), 'ContentType', 'vector');
cd(parentDir)
%% 2. amplification for 4 genes (supplementary figure), page width
xplotNames=            {'UBC\namplifier set 9\nAtto488 readout',...
    'ITGA3\namplifier set 10\nCy3 readout',...
    'FN1\namplifier set 5\nAlexa594 readout',...
    'MITF\namplifier set 12\nAtto647N readout'};
plotChannelBaseNames={'gfp','tmr' ,'alexa' ,'cy' };
plotPositionNums=     [ 1      2     3     4      ];
numPlots=length(plotChannelBaseNames);

figureWidthPoints=540;
figureHeightPoints=180;

% condID to plot
plotCondIDs=[1 2 3 4 5 6];

assert(all(strcmp(Tcond.SpecialClick(ismember(Tcond.condID,plotCondIDs)),''))) %check for consistency

cd(parentDir)
figure(2);
fh=gcf; clf; hold on;
set(fh,'Position',[1500 400 figureWidthPoints figureHeightPoints],'Units','points')
for iPlot=1:numPlots
    channelBaseName=plotChannelBaseNames{iPlot};
    plotName=plotNames{iPlot};
    positionNum=plotPositionNums(iPlot);
    ampRoundsOrdered=Tcond.AmpRound(ismember(Tcond.condID,plotCondIDs));
    
    indTspots=all([strcmp(Tspots.channelBaseName,channelBaseName),ismember(Tspots.condID,plotCondIDs)],2);
    X=Tspots.AmpRound(indTspots);
    %Xcat=Tspots.condName(indTspots);
    %Xcat=replace(Xcat,'R','');
    %Xcat=replace(Xcat,{'R','_NoCu'},'');
    Y=Tspots.intensitiesNorm(indTspots);
    
    % create subplots within the figure, 1x4 arrangement, and select one
    subplot(1,4,positionNum)
    
    ax(iPlot)=gca;
    ax(iPlot).PositionConstraint = 'outerposition';
    % set position of current axis
    XaxesWidth=0.225;
    XleftSpace=0.035;
    XbetweenAxesSpace=0.01;
    YbotSpace=0.12;
    YaxesHeight=0.87;
    XPosAxes= XleftSpace + (iPlot-1)*(XaxesWidth + XbetweenAxesSpace);
    ax(iPlot).Position=[XPosAxes YbotSpace XaxesWidth YaxesHeight];
    
    % plot violins
    vs=violinplot(Y,X,'ViolinColor',clickPositiveColor,'ShowData',false,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',BoxWidth,'BoxColor',BoxColor);
    
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs)
            vs(iViolin).WhiskerPlot.LineStyle='none';
        end
    end
    
    % median marker
    for iViolin=1:length(vs)
        vs(iViolin).MedianPlot.Marker=MedianMarker;
        vs(iViolin).MedianPlot.LineWidth=MedianLineWidth;
        vs(iViolin).MedianPlot.SizeData=MedianSizeData;
    end
    
    % log scale
    ax(iPlot).YScale='log';
    
    % XLim and YLim
    ax(iPlot).XLim=[0 11.5];
    ax(iPlot).YLim=[0.5 220]; % set equal to reference
    
    % make violin X's be centered around AmpRound (1,2,4,6,8,10)
    XDistanceToAddArray=nan(1,length(vs));
    for iViolin=1:length(vs)
        thisAmpRound=ampRoundsOrdered(iViolin);
        XDistanceToAddArray(iViolin)=thisAmpRound - iViolin;
    end
    vs=XshiftViolinplot(vs,XDistanceToAddArray);
    
    % X Tick Labels
    ax(iPlot).XAxis.TickValues=ampRoundsOrdered;
    
    % Y Tick Labels (turn off for right-side subplots)
    if any(positionNum==[2 3 4]')
        ax(iPlot).YTickLabel=[];
    end
    
    
    
    % Y Ticks
    ax(iPlot).YTick=[1 2 4 8 16 32 64 128];
    ax(iPlot).YMinorTick='off';
    
    % Font size
    ax(iPlot).FontSize=fontSize;
    
    % make axis background invisible 
    set(ax(iPlot), 'color', 'none');

    
    % title
    text(.02, 0.9, sprintf(plotName), 'Units', 'normalized', 'FontSize',titleFontSize,'FontWeight','normal')
    
    % medians for regression
    Ymedians=nan(length(plotCondIDs),1);
    for iCondID=1:length(plotCondIDs)
        condID=plotCondIDs(iCondID);
        YdataThisGeneThisCond=Tspots.intensitiesNorm(all([strcmp(Tspots.channelBaseName,channelBaseName),Tspots.condID==condID],2));
        Ymedians(iCondID)=median(YdataThisGeneThisCond);
    end
    if plot_medians_regression_line
        XdataReg=ampRoundsOrdered(2:end);
        YdataReg=Ymedians(2:end);
        %f = fit(XdataReg,YdataReg,'power2');
        coeff=polyfit(XdataReg',log2(YdataReg'),1); % use data before the jitter
        m=coeff(1); b=coeff(2);
        %XregLim=[ampRoundsOrdered(2) ampRoundsOrdered(end)];
        XregLim=XdataReg;
        YregLim=2.^(m*XregLim+b);
        plot(ax(iPlot),XregLim,YregLim,'-','Color','k')
        regressionText=sprintf('y=%.3f\\times%.3f^{Round}',2^b,2^m);
        %text(ax(iPlot),.30,0.17,regressionText,'FontSize',regressionFontSize,'Interpreter','none','Units', 'normalized')
        text(ax(iPlot),.37,0.17,regressionText,'FontSize',regressionFontSize,'Units', 'normalized')
    end
    
    % medians text labels
    if medians_textlabels
        for iCondID=1:length(plotCondIDs)
            ampRound=ampRoundsOrdered(iCondID);
            thisMedian=Ymedians(iCondID);
            thisMedianLabelOffsetX=zeros(size(thisMedian));
            thisMedianLabelOffsetY=zeros(size(thisMedian));
            if and(iPlot==4,iCondID==1)
                thisMedianLabelOffsetX(1)=-0.25;
                thisMedianLabelOffsetY(1)=-0.3;
            end
            text(ax(iPlot),ampRound+0.4+thisMedianLabelOffsetX,thisMedian+thisMedianLabelOffsetY,num2str(thisMedian,'%0.1f'),'FontSize',regressionFontSize)
        end
    end
    
    % xlabel for whole figure
    if positionNum==2
        xlabelText='Round of amplification';
        text(0.7, -0.1, xlabelText, 'Units', 'normalized', 'FontSize',axisLabelFontSize,'FontWeight','normal')
    end
    
    % ylabel for whole figure
    if positionNum==1
        ylabelText= text(-0.13,0.2,'Spot intensity, normalized','Units','normalized','FontSize',axisLabelFontSize,'FontWeight','normal');
        set(ylabelText,'Rotation',90)
    end
    
end

%fh=packSubplotsTogether(fh,1.25,1.25);

cd(parentDirForExport)
if ~isfolder(plotDir)
    mkdir(plotDir)
end
exportgraphics(gcf, fullfile(plotDir,'E159_amplification_4genes.eps'), 'ContentType', 'vector');
cd(parentDir)


%% 3. +/-click for 4 genes (supplementary figure), page width
plotNames=            {'UBC\namplifier set 9\nAtto488 readout',...
    'ITGA3\namplifier set 10\nCy3 readout',...
    'FN1\namplifier set 5\nAlexa594 readout',...
    'MITF\namplifier set 12\nAtto647N readout'};
plotChannelBaseNames={'gfp','tmr' ,'alexa' ,'cy' };
plotPositionNums=     [ 1      2     3     4      ];
numPlots=length(plotChannelBaseNames);

figureWidthPoints=540;
figureHeightPoints=180;


plotCondIDs=[2 7 5 8];
ampRoundsOrdered=Tcond.AmpRound(ismember(Tcond.condID,plotCondIDs));
specialClicks={'','NoCu','','NoCu'}';

cd(parentDir)
figure(3); clf; hold on;
fh=gcf;
set(fh,'Position',[1500 100 figureWidthPoints figureHeightPoints],'Units','points')

clear ax
for iPlot=1:numPlots
    channelBaseName=plotChannelBaseNames{iPlot};
    plotName=plotNames{iPlot};
    positionNum=plotPositionNums(iPlot);
    
    indTspots=all([strcmp(Tspots.channelBaseName,channelBaseName),ismember(Tspots.condID,plotCondIDs)],2);
    Xcat=Tspots.condName(indTspots);
    assert(isequal(unique(Xcat),{'R2';'R2_NoCu';'R8';'R8_NoCu'}));
    Y=Tspots.intensitiesNorm(indTspots);
    
    subplot(1,4,positionNum)
    
    % set axes size
    ax(iPlot)=gca;
    ax(iPlot).PositionConstraint = 'outerposition';
    XplotWidth=0.225;
    XleftSpace=0.035;
    XbetweenSpace=0.01;
    YbotSpace=0.12;
    YplotHeight=0.87;
    XPosAxes= XleftSpace + (iPlot-1)*(XplotWidth + XbetweenSpace);
    ax(iPlot).Position=[XPosAxes YbotSpace XplotWidth YplotHeight];
    
    ax(iPlot).YScale='log';
    
    ax(iPlot).XLim=[0 5];
    
    vs=violinplot(Y,Xcat,'ViolinColor',[0 0 1],'ShowData',false,'MedianColor',MedianColor,'ViolinAlpha',ViolinAlpha,'BoxWidth',0.02,'BoxColor',[0.4 0.4 0.4]);
    
    % XLim and YLim
    ax(iPlot).XLim=[0 5];
    ax(iPlot).YLim=[0.5 220]; % set equal to reference
    
    
    
    % recolor the violins
    for iViolin=1:length(vs)
        specialClick=specialClicks{iViolin};
        switch specialClick
            case ''
                vs(iViolin).ViolinPlot.FaceColor=clickPositiveColor;
            case 'NoCu'
                vs(iViolin).ViolinPlot.FaceColor=clickNegativeColor;
            otherwise
                error('case not handled')
        end
    end
    
    % delete the whiskers
    if removeWhiskers
        for iViolin=1:length(vs)
            vs(iViolin).WhiskerPlot.LineStyle='none';
        end
    end
    
    
    % median marker
    for iViolin=1:length(vs)
        vs(iViolin).MedianPlot.Marker=MedianMarker;
        vs(iViolin).MedianPlot.LineWidth=MedianLineWidth;
        vs(iViolin).MedianPlot.SizeData=MedianSizeData;
    end
    
    % shift violins' X's
    vs=XshiftViolinplot(vs,[0.15 -0.15 0.15 -0.15]);
    
    % X Ticks
    ax(iPlot).XAxis.TickLabelInterpreter='none';
    ax(iPlot).XTick=[1.5 3.5];
    ax(iPlot).XTickLabel=[2 8];
    if any(positionNum==[]') %put iPlots that shouldn't have x tick labels
        ax(iPlot).XTickLabel=[];
    end
    
    % Y Tick Labels (turn off for right-side subplots)
    if any(positionNum==[2 3 4]')
        ax(iPlot).YTickLabel=[];
    end
    ax(iPlot).YTick=[1 2 4 8 16 32 64 128];
    ax(iPlot).YMinorTick='off';
    
    % font size
    ax(iPlot).FontSize=fontSize;
    
        % make axis background invisible 
    set(ax(iPlot), 'color', 'none');

    
    text(.02, 0.9, sprintf(plotName), 'Units', 'normalized', 'FontSize',titleFontSize,'FontWeight','normal')
    
    % xlabel for whole figure
    xlabelText='Round of amplification';
    if iPlot==2
        text(0.7, -0.1, xlabelText, 'Units', 'normalized', 'FontSize',axisLabelFontSize,'FontWeight','normal')
    end
    
    % ylabel for whole figure
    if iPlot==1
        htext = text(-0.13,0.2,'Spot intensity, normalized','Units','normalized','FontSize',axisLabelFontSize,'FontWeight','normal');
        set(htext,'Rotation',90)
    end
    
    
    % medians text labels
    if medians_textlabels
        
        % medians for regression
        Ymedians=nan(length(plotCondIDs),1);
        for iCondID=1:length(plotCondIDs)
            condID=plotCondIDs(iCondID);
            YdataThisGeneThisCond=Tspots.intensitiesNorm(all([strcmp(Tspots.channelBaseName,channelBaseName),Tspots.condID==condID],2));
            Ymedians(iCondID)=median(YdataThisGeneThisCond);
        end
        for iCondID=1:length(plotCondIDs)
            ampRound=ampRoundsOrdered(iCondID);
            thisMedian=Ymedians(iCondID);
            text(ax(iPlot),iCondID+0.3,thisMedian*1.3,num2str(thisMedian,'%0.1f'),'FontSize',regressionFontSize)
        end
    end
    
    % legend
    if iPlot==4
        text(0.7, 0.85,'(+) Click','Color',clickPositiveColor,'FontSize',fontSize,'Units','normalized')
        text(0.7, 0.78,'(-) Click','Color',clickNegativeColor,'FontSize',fontSize,'Units','normalized')
    end
    
end

cd(parentDirForExport)
%saveas(gcf,fullfile(plotDir,filesep,'E159_click.svg'))
exportgraphics(gcf,fullfile(plotDir,filesep,'E159_click_4genes.eps'), 'ContentType', 'vector')
cd(parentDir)











