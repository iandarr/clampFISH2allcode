% display zoom-in of spots to show that they don't grow in physical size

parentDir='/Volumes/IAND_08'; % for raw data
Tcond=readtable('/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/experiments/E159part1rep2_ampclick_Tcond.xlsx');
parentDirForExtractedData='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
extractedDataDirAllConditions='paper/extractedData/E159part1rep2_ampclick_100X';
Tspots=readtable(fullfile(parentDirForExtractedData,filesep,extractedDataDirAllConditions,'Tspots1All.csv'));


plotDir=[parentDirForExtractedData,filesep,'paper/plots/E159part1rep2_E166_spotSize/']; % don't need to export to exampleImages as individual tiffs since these are very low-res.
plotName='example_spot_ZoomIns';


condIDsToShow=[1:6]';
numPlotCols=length(condIDsToShow);
numSpotsToExractPerCond=3;

%channel='cyb';
channel='gfpb';
halfWidthOfRegionToExtract=3; % pixels. cropped image will have a width of 1+2*halfWidthOfRegion
assert(rem(halfWidthOfRegionToExtract,1)==0)
rng(0)
croppedSpotImgs=zeros([1+halfWidthOfRegionToExtract*2,1+halfWidthOfRegionToExtract*2,numSpotsToExractPerCond,numPlotCols],'uint16');
TspotsExtracted=table();

for i=1:numPlotCols
    colNum=i;
    condID=condIDsToShow(i);
    rowTcond=find(Tcond.condID==condID);
    AmpRound=Tcond.AmpRound(rowTcond);
    relpathRawDataDir=Tcond.Tiff_100X_dir{rowTcond};

    TspotsThisCond=Tspots(all([strcmp(Tspots.channel,channel), Tspots.condID==condID],2),:);
    %indSpots=randperm(height(TspotsThisCond),numSpotsToExractPerCond);
    % get median intensity spots
    [~,indSorted]=sort(TspotsThisCond.amplitudeFitted);
    midIndex=round(length(indSorted)/2);
    midIndices=midIndex:(midIndex+numSpotsToExractPerCond-1);
    indSpots=indSorted(midIndices);
    TspotsExctractedThisCond=TspotsThisCond(indSpots,:);
    arrayNums=TspotsExctractedThisCond.arrayNum;
    objNums=TspotsExctractedThisCond.objNum;
    XYs=[TspotsExctractedThisCond.X,TspotsExctractedThisCond.Y];
    Zs=TspotsExctractedThisCond.Z;
    %fullpathRawDataDir=
    TspotsExctractedThisCond.plotColNum=repmat(i,numSpotsToExractPerCond,1);
    TspotsExctractedThisCond.condID=repmat(condID,numSpotsToExractPerCond,1);
    TspotsExctractedThisCond.spotInd=[1:numSpotsToExractPerCond]';
    TspotsExtracted=[TspotsExtracted;TspotsExctractedThisCond];
    % get cropped images
    for iSpot=1:numSpotsToExractPerCond
        arrayNum=arrayNums(iSpot);
        objNum=objNums(iSpot);
        zplane=Zs(iSpot);
        fileSuffixStr=sprintf('%03d',arrayNum);
        fullpathDir=fullfile(parentDir,filesep,relpathRawDataDir,filesep);
        fullpathFilename=fullfile(parentDir,filesep,relpathRawDataDir,filesep,[channel,fileSuffixStr,'.tif']);
        ColRange=[XYs(iSpot,1)-halfWidthOfRegionToExtract, XYs(iSpot,1)+halfWidthOfRegionToExtract];
        RowRange=[XYs(iSpot,2)-halfWidthOfRegionToExtract, XYs(iSpot,2)+halfWidthOfRegionToExtract];

        getSpots=false; % don't need to have this be true. Already have Tspots.
        threshold=0; % just choose something
        [imgFISHSource,croppedMask,spotRowColZ_check]=getRajlabimagetoolsCroppedCellImage(fullpathDir,arrayNum,objNum,channel,threshold,zplane,getSpots);
        croppedSpotImgs(:,:,iSpot,i)=imgFISHSource(RowRange(1):RowRange(2),ColRange(1):ColRange(2));
        %croppedSpotImgs(:,:,iSpot)=imread(fullpathFilename,zplane,'PixelRegion',{RowRange,ColRange});
    end
end

%% show croppped spot images
%TspotsExtracted.spotInd=repmat([1:numSpotsToExractPerCond]',numPlotCols,1);
spotsIndToPlot=3;
TspotsPlotted=TspotsExtracted(TspotsExtracted.spotInd==spotsIndToPlot,:);

NmPerPx=65;

%numSpotsToShowPerCond=length(spotsIndToPlotPerCond);
%numPlotCols=length(condIDsToShow);  
assert(all(ismember(spotsIndToPlot,1:numSpotsToExractPerCond)))

f=figure(1);clf;
FigurePosition=[50 50 300 150];
f.Position=FigurePosition;

%t=tiledlayout(numSpotsToShowPerCond,numPlotCols);
t=tiledlayout(2,numPlotCols);
t.TileIndexing='rowmajor';
t.TileSpacing='compact';
t.Padding='normal';
%xlabel(t,'Round of amplification','FontSize',8)
%ylabel(t,sprintf('standard deviation\nof gaussian fit of spots (nm)'))
t.YLabel.FontSize=8;
t.XLabel.Color=[0 0 0]; t.YLabel.Color=[0 0 0];


minIntensityAllImgs=min(croppedSpotImgs(:,:,spotsIndToPlot,:),[],'all');
maxIntensityAllImgs=max(croppedSpotImgs(:,:,spotsIndToPlot,:),[],'all');
normAmplitude=TspotsPlotted.amplitudeFitted(TspotsPlotted.AmpRound==1);

for i=1:numPlotCols
    condID=condIDsToShow(i);
    rowTcond=find(Tcond.condID==condID);
    AmpRound=Tcond.AmpRound(rowTcond);

    colNum=i;
    img=croppedSpotImgs(:,:,spotsIndToPlot,i);
    CLimLow=min(croppedSpotImgs(:,:,spotsIndToPlot,i),[],'all');
    CLimHigh=max(croppedSpotImgs(:,:,spotsIndToPlot,i),[],'all');

    % FIRST ROW: same contrasting
    rowNum=1;
    axisNum=numPlotCols*(rowNum-1)+colNum;
    ax=nexttile(axisNum);
    imagesc(ax,img,[minIntensityAllImgs,maxIntensityAllImgs])
    axis equal; axis square;
    ax.XAxis.TickLabels=[];ax.YAxis.TickLabels=[];
    colormap gray;
    textY=halfWidthOfRegionToExtract*2+2.1;
    text(0.5,textY,sprintf(' %i-%i',minIntensityAllImgs,maxIntensityAllImgs),'Color',0.3*[1 1 1],'FontSize',6) % contrast
    ax.XLim=[0.5 halfWidthOfRegionToExtract*2+1.5];
    ax.YLim=[0.5 halfWidthOfRegionToExtract*2+1.5];
    ax.XTick=[];ax.YTick=[];
    %set(ax,'XColor','none');set(ax,'YColor','none');
    %ax.XAxis.Color='none';ax.YAxis.Color='none'
    %box off;
    if i==1
        %ylabel(sprintf('equal\ncontrasting'),'FontSize',8)
        ylabel(sprintf('contrasted\nequally'),'FontSize',8,'Color',[0 0 0])
    end
    
    % spot info
    amplitudeFittedNormalized=TspotsPlotted.amplitudeFitted(TspotsPlotted.condID==condID)/normAmplitude;
    sigmaFittedNm=TspotsPlotted.sigmaFitted(TspotsPlotted.condID==condID)*NmPerPx;
    %title(sprintf('round %i\n%.2f\n%.0fnm',AmpRound,amplitudeFittedNormalized,sigmaFittedNm),'FontSize',8,'FontWeight','normal')
    %title(sprintf('round %i',AmpRound),'FontSize',8,'FontWeight','normal','Color',[0 0 0])
    text(0.5,0,sprintf('  round %i',AmpRound),'FontSize',8,'FontWeight','normal','Color',[0 0 0])

    % SECOND ROW: column-specific contrasting
    rowNum=2;
    axisNum=numPlotCols*(rowNum-1)+colNum;
    ax=nexttile(axisNum);
    imagesc(ax,img,[CLimLow,CLimHigh])
    axis equal; axis square;
    ax.XAxis.TickLabels=[];ax.YAxis.TickLabels=[];
    colormap gray;
    textY=halfWidthOfRegionToExtract*2+2.1;
    text(0.5,textY,sprintf('  %i-%i',CLimLow,CLimHigh),'Color',0.3*[1 1 1],'FontSize',6) % contrast
    ax.XLim=[0.5 halfWidthOfRegionToExtract*2+1.5];
    ax.YLim=[0.5 halfWidthOfRegionToExtract*2+1.5];
    ax.XTick=[];ax.YTick=[];
    %set(ax,'XColor','none');set(ax,'YColor','none');
    %ax.XAxis.Color='none';ax.YAxis.Color='none'
    %box off;
    if i==1
        %ylabel(sprintf('individual\ncontrasting'),'FontSize',8,'Color',[0 0 0])
        ylabel(sprintf('contrasted\nindividually'),'FontSize',8,'Color',[0 0 0])
    end
    %end

    % info about spots
    textLineSpacing=1.5;
    text(0.5,textY+textLineSpacing*1,sprintf('   %.2f',amplitudeFittedNormalized),'Color',[0 0 0],'FontSize',8) % contrast
    text(0.5,textY+textLineSpacing*2,sprintf('  %.0fnm',sigmaFittedNm),'Color',[0 0 0],'FontSize',8) % contrast
end

exportgraphics(f,fullfile(plotDir,filesep,[plotName,'.eps']), 'ContentType', 'vector')
