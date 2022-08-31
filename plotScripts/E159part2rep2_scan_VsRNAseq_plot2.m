% E159part2_scan_VsRNAseq_plot2
% linear plot with zoom-in
% high-throughput scan, average FISH vs. bulk RNAseq
% for:
%   - drug-naive WM989 A6-G3 (wells A1,A2,A3,B1,B2)
%   - 1uM vemurafenib-resistant WM989 A6-G3 RC4 (well B3)

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
tableDir='paper/extractedData/E159part2rep2_scan/VsRNAseq/';
plotOutDir='paper/plots/E159part2rep2_scan/VsRNAseq/';


T_RNAseq_naive=readtable(fullfile(parentDir,tableDir,'tpm_RNAseq_naiveSubset.csv'))
T_RNAseq_resistant=readtable(fullfile(parentDir,tableDir,'tpm_RNAseq_resistantSubset.csv'))
T_clamp2_naive=readtable(fullfile(parentDir,tableDir,'clamp2_meanCount_naive.csv'))
T_clamp2_resistant=readtable(fullfile(parentDir,tableDir,'clamp2_meanCount_resistant.csv'))

%%
plotOutName='mean_clamp2_vs_RNAseq_linearRegression';
plot_regression=true;

plot_inset=true;

fontSize=6;
fontSizeSmaller=6;
plot_log2_values=true;
plot_regression_line=true;

markerSize=6;
MarkerFaceColor=[0 109 219]/255;
MarkerFaceAlpha=0.6;
% check variable ordering
isequal(T_RNAseq_naive.Properties.VariableNames,T_RNAseq_resistant.Properties.VariableNames,T_clamp2_naive.Properties.VariableNames,T_clamp2_resistant.Properties.VariableNames)

geneSymbols=T_RNAseq_naive.Properties.VariableNames;
assert(length(geneSymbols)==10)

XY_naive=    [T_clamp2_naive{:,:}' ,T_RNAseq_naive{:,:}'];
XY_resistant=[T_clamp2_resistant{:,:}' ,T_RNAseq_resistant{:,:}'];

f=figure(1); clf
set(f,'Position',[1000 200 275 200],'Units','points')
t=tiledlayout(1,2);
t.TileSpacing='compact';
t.Padding='compact';
for iPlot=1:2
    ax=nexttile(t,iPlot); hold on;

    if iPlot==1
        XY=XY_naive;
        cellType={'drug-naive','WM989 A6-G3'};
        ax.XLim=[-5 80];
        ax.YLim=[-40 650];
    elseif iPlot==2
        XY=XY_resistant;
        cellType={'drug-resistant','WM989 A6-G3 RC4'};
        ax.XLim=[-50 800];
        ax.YLim=[-100 2200];
    else
        error('')
    end

    scatter(ax,XY(:,1),XY(:,2),markerSize,'filled','MarkerFaceColor',MarkerFaceColor,'MarkerFaceAlpha',MarkerFaceAlpha);

    ax.FontSize=fontSize;

    %ax.XLim=[10^-1 10^4];
    %ax.YLim=[10^-1 10^4];


    title(cellType,'FontWeight','normal','FontSize',fontSize)

    if plot_regression_line
        hold on
        mdl=fitlm(XY(:,1),XY(:,2));
        %mdl=fitlm(X,Y);
        Tcoeff=mdl.Coefficients;
        m=Tcoeff.Estimate('x1'); % slope
        b=Tcoeff.Estimate('(Intercept)'); % intercept
        Rsq=mdl.Rsquared.Ordinary;
    end

    % plot regression
    if plot_regression
        hreg=fplot(ax,@(x) m*x + b);
        hreg.Color='k';
    end

    % regressionText
    if sign(b)==1
        AddOrSubtractB='+';
    elseif sign(b)==-1
        AddOrSubtractB='-';
    else
        error('huh')
    end

    regressionText=sprintf('R^2=%.3f\ny=%.2fx %s %.1f',Rsq,m,AddOrSubtractB,abs(b));
    text(0.6*ax.XLim(2),0.02*ax.YLim(2),regressionText,'FontSize',fontSizeSmaller)
    regressionText
    % data labels

    %                       above regression line           below regression line
    XlabelOffset=-1  *(XY(:,2)>=m*XY(:,1)+b)    +  0.25*(XY(:,2)<m*XY(:,1)+b);
    YlabelOffset= .5  *(XY(:,2)>=m*XY(:,1)+b)    + -0.25*(XY(:,2)<m*XY(:,1)+b);

    geneSymbolsItalicized=join([repmat({'\it'},1,10);geneSymbols],'',1);
    


    % plot inset
    if plot_inset
     if iPlot==1
         axIn=axes('Position',[.10 .6 .21 .21]);
         xFact=10;
         yFact=30;
         %axIn=axes('Position',[.3 .2 .15 .15]);
         elseif iPlot==2
             axIn=axes('Position',[.6 .6 .21 .21]);
         xFact=6;
         yFact=13;
         else
        error('')
     end

        axIn.XLim=[ax.XLim(1)/(0.5*xFact) ax.XLim(2)/xFact];
        axIn.YLim=[ax.YLim(1)/(0.5*yFact) ax.YLim(2)/yFact];
      idx=XY(:,1)<=axIn.XLim(2) & XY(:,2)<=axIn.YLim(2);
      sum(idx)
           scatter(axIn,XY(idx,1),XY(idx,2),markerSize,'filled','MarkerFaceColor',MarkerFaceColor,'MarkerFaceAlpha',MarkerFaceAlpha);
        axIn.XLim=[ax.XLim(1)/(0.5*xFact) ax.XLim(2)/xFact];
        axIn.YLim=[ax.YLim(1)/(0.5*yFact) ax.YLim(2)/yFact];
        text(XlabelOffset(idx)+XY(idx,1),YlabelOffset(idx)+XY(idx,2),geneSymbolsItalicized(idx),'FontSize',fontSizeSmaller)
    axIn.FontSize=fontSizeSmaller;
    box on
        % non-inset labels
        text(ax,XlabelOffset(~idx)+XY(~idx,1),YlabelOffset(~idx)+XY(~idx,2),geneSymbolsItalicized(~idx),'FontSize',fontSizeSmaller,'Color','k')

    else
        % all labels
        text(ax,XlabelOffset+XY(:,1),YlabelOffset+XY(:,2),geneSymbolsItalicized,'FontSize',fontSizeSmaller,'Color','k')
    end
    
    axIn.XAxis.Color=[0 0 0];
    axIn.XAxis.Label.Color=[0 0 0];
    ax.XAxis.Color=[0 0 0];
    ax.XAxis.Label.Color=[0 0 0];
    
    axIn.YAxis.Color=[0 0 0];
    axIn.YAxis.Label.Color=[0 0 0];
    ax.YAxis.Color=[0 0 0];
    ax.YAxis.Label.Color=[0 0 0];
end



%if plot_log2_values
%    xlabel(t,'bulk RNA-seq log_2(transcripts per million, tpm)','FontSize',fontSize)
%    ylabel(t,'mean clampFISH 2.0 log_2(RNA spots)','FontSize',fontSize)
%else
xlabel(t,'mean clampFISH 2.0 RNA spots','FontSize',fontSize,'Color','k')
ylabel(t,'bulk RNA-seq transcripts per million (tpm)','FontSize',fontSize,'Color','k')
%end
exportgraphics(f, fullfile(plotOutDir,filesep, [plotOutName,'.eps']), 'ContentType', 'vector'); % with spots
