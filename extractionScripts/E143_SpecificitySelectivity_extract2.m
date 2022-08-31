% E143_ampScreen_Colocalization
%% E143_ampScreen_plot
% amplifier screen

thresholdList=... %This threshold list takes about 13hr to complete
    [ 100  100  100  100  100  100  100  100  100   125   250  200  100   250  150,...
    450  150  400  250  250  300  200  200  400   500  1000  800  160  1000  650,... %GFP
    1900  270 1200  700  650 1000  850 1200 2000  3100  2900 1800  300  2300 1700,...%EGFR
    NaN  NaN  NaN]';

condID2threshForExtraction=containers.Map(1:45,thresholdList(1:45));
%% Go to parent folder
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
% Read in data
Tcond=readtable(fullfile(parentDir,'paper/experiments/E143_amplifierScreen_Tcond.xlsx'));

extractedDataDir='paper/extractedData/E143_screen/60X';

Tspots=readtable(fullfile(parentDir,filesep,extractedDataDir,filesep,'TspotsAll.csv'));
Tspots.isGood=logical(Tspots.isGood);
%Tcell=readtable(fullfile(parentDir,filesep,extractedDataDir,filesep,'Tcell2.csv'));
condIDList=[1:30]';

% export this new extractedData to clamp2rev
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2rev';
extractedDataOutputDir='paper/extractedData/E143_SensitivitySelectivity/';
%% True positive and false false positive clampFISH 2.0 spots
smFISHchannel='tmr';
clamp2channel='cya';
colocDistXY=2.8; % xy pixels. 60X imaging at 2x2 binning
colocDistZ=1;

xShift=10;
yShift=0;

numCond=length(condIDList);
numThresholds=5;
threshFoldChange=16;
Tres=array2table(nan(numCond*numThresholds,15),'VariableNames',{'isShifted','condID','condName','primary','amplifierSeries','threshold','threshExtracted','propClampColoc','propSmFISHColoc','numClampColoc','numSmFISHColoc','numClampTot','numSmFISHTot','xShift','yShift'});
Tres=convertvars(Tres,{'condName','primary'},'string');
Tres=convertvars(Tres,{'condName','primary'},'char');
TresSHIFTED=Tres;
iRes=0;
for i=1:numCond
    condID=condIDList(i);
    condName=Tcond.condName{Tcond.condID==condID};
    primary=Tcond.primary{Tcond.condID==condID};
    amplifierSeries=Tcond.amplifierSeries(Tcond.condID==condID);

    fprintf('starting colocalization analysis for condID=%i (%s) at %s\n',condID,condName,datetime)
    minThresh=condID2threshForExtraction(condID);
    highThresh=threshFoldChange*minThresh;

    thresholds=round(logspace(log10(minThresh), log10(highThresh), numThresholds));
    thresholds(2:end-1)=round(thresholds(2:end-1),-2);
    %LowHighThresh=[400 2000];

    % CONDITION-SPECIFIC
    %minThresh=LowHighThresh(1);
    TspotsClamp=Tspots(all([Tspots.condID==condID,Tspots.isGood,strcmp(Tspots.channel,clamp2channel),Tspots.intensities>=minThresh],2),:);
    TspotsSmFISH=Tspots(all([Tspots.condID==condID,Tspots.isGood,strcmp(Tspots.channel,smFISHchannel)],2),:);
    % CELL-SPECIFIC
    C=struct();
    cellIDlist=unique(TspotsClamp.cellID);
    assert(isequal(cellIDlist,unique(TspotsSmFISH.cellID)))
    numCells=length(cellIDlist);
    for iC=1:numCells
        cellID=cellIDlist(iC);
        C(iC).cellID=cellID;


        idxClampThisCell=TspotsClamp.cellID==cellID;
        idxSmFISHThisCell=TspotsSmFISH.cellID==cellID;
        C(iC).numSmFISHSpots=sum(idxSmFISHThisCell);
        % inefficient to do it in one big matrix, but also easy:

        C(iC).intensitiesClamp=TspotsClamp{idxClampThisCell,{'intensities'}}; % col vect
        xyzClamp=TspotsClamp{idxClampThisCell,{'X','Y','Z'}}; % col vect
        xyzSmFISH=TspotsSmFISH{idxSmFISHThisCell,{'X','Y','Z'}};

        distxyM=sqrt((xyzClamp(:,1) - xyzSmFISH(:,1)').^2 +...
            (xyzClamp(:,2) - xyzSmFISH(:,2)').^2); % one row per clamp2 spot, one column per smFISH spot
        distzM=sqrt((xyzClamp(:,3) - xyzSmFISH(:,3)').^2);
        C(iC).isColocM=(distxyM<=colocDistXY) & (distzM<=colocDistZ);

        % pixel shifted data
        distxyMSHIFTED=sqrt(((xyzClamp(:,1) +xShift) - xyzSmFISH(:,1)').^2 +...
            ((xyzClamp(:,2) +yShift) - xyzSmFISH(:,2)').^2); % one row per clamp2 spot, one column per smFISH spot
        C(iC).isColocMSHIFTED=(distxyMSHIFTED<=colocDistXY) & (distzM<=colocDistZ);
    end

    numSmFISHSpotsAllCells=height(TspotsSmFISH);
    for ii=1:numThresholds
        threshold=thresholds(ii);
        
        numClampSpotsAllCells=0;
        numColocClampAllCells=0;
        numColocSmFISHAllCells=0;
        
        numColocClampAllCellsSHIFTED=0;
        numColocSmFISHAllCellsSHIFTED=0;
        for iC=1:numCells
            isAboveThreshM=repmat(C(iC).intensitiesClamp>=threshold,1,C(iC).numSmFISHSpots);
            numClampSpots=sum(C(iC).intensitiesClamp>=threshold);

            % real data
            isColocThisThreshM=C(iC).isColocM & isAboveThreshM;
            idxClampColocWithSmFISH=any(isColocThisThreshM,2); % col vect
            idxSmFISHColocWithClamp=any(isColocThisThreshM,1); % row vect
            numClampSpotsAllCells=numClampSpotsAllCells+numClampSpots;
            numColocClampAllCells=numColocClampAllCells + sum(idxClampColocWithSmFISH);
            numColocSmFISHAllCells=numColocSmFISHAllCells+ sum(idxSmFISHColocWithClamp);

            % pixel-shifted data
            isColocThisThreshMSHIFTED=C(iC).isColocMSHIFTED & isAboveThreshM;
            idxClampColocWithSmFISHSHIFTED=any(isColocThisThreshMSHIFTED,2); % col vect
            idxSmFISHColocWithClampSHIFTED=any(isColocThisThreshMSHIFTED,1); % row vect
            
            numColocClampAllCellsSHIFTED=numColocClampAllCellsSHIFTED+sum(idxClampColocWithSmFISHSHIFTED);
            numColocSmFISHAllCellsSHIFTED=numColocSmFISHAllCellsSHIFTED+sum(idxSmFISHColocWithClampSHIFTED);
        end
        propClampColoc=numColocClampAllCells/numClampSpotsAllCells;
        propSmFISHColoc=numColocSmFISHAllCells/numSmFISHSpotsAllCells;

        propClampColocSHIFTED=numColocClampAllCellsSHIFTED/numClampSpotsAllCells;
        propSmFISHColocSHIFTED=numColocSmFISHAllCellsSHIFTED/numSmFISHSpotsAllCells;

        % fill in results table
        iRes=iRes+1;
        Tres.xShift(iRes)=0;
        Tres.yShift(iRes)=0;
        Tres.condID(iRes)=condID;
        Tres.condName(iRes)={condName};
        Tres.primary(iRes)={primary};
        Tres.amplifierSeries(iRes)=amplifierSeries;
        Tres.threshold(iRes)=threshold;
        Tres.threshExtracted(iRes)=minThresh;
        Tres.propClampColoc(iRes)=propClampColoc;
        Tres.propSmFISHColoc(iRes)=propSmFISHColoc;
        Tres.numClampColoc(iRes)=numColocClampAllCells;
        Tres.numSmFISHColoc(iRes)=numColocSmFISHAllCells;
        Tres.numClampTot(iRes)=numClampSpotsAllCells;
        Tres.numSmFISHTot(iRes)=numSmFISHSpotsAllCells;
        Tres.isShifted(iRes)=false;
        % fill in shifted results table
        TresSHIFTED(iRes,{'condID','condName','primary','amplifierSeries','threshold','threshExtracted','numClampTot','numSmFISHTot'})=Tres(iRes,{'condID','condName','primary','amplifierSeries','threshold','threshExtracted','numClampTot','numSmFISHTot'});
        TresSHIFTED.xShift(iRes)=xShift;
        TresSHIFTED.yShift(iRes)=yShift;
        TresSHIFTED.propClampColoc(iRes)=propClampColocSHIFTED;
        TresSHIFTED.propSmFISHColoc(iRes)=propSmFISHColocSHIFTED;
        TresSHIFTED.numClampColoc(iRes)=numColocClampAllCellsSHIFTED;
        TresSHIFTED.numSmFISHColoc(iRes)=numColocSmFISHAllCellsSHIFTED;
        TresSHIFTED.numClampTot(iRes)=numClampSpotsAllCells;
        TresSHIFTED.numSmFISHTot(iRes)=numSmFISHSpotsAllCells;
        TresSHIFTED.isShifted(iRes)=true;
    end
end

Tres
TresSHIFTED

TresOut=[Tres;TresSHIFTED];
writetable(TresOut,fullfile(parentDir,filesep,extractedDataOutputDir,['TsensitivitySelectivity.csv']))