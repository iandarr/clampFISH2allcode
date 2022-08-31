% E143_ampScreen_extract
% amplifier screen

%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E143_amplifierScreen_Tcond.xlsx');



%% Segmentation
%navigate to a subfolder and run this. Segment cells.

%improc2.segmentGUI.SegmentGUI

%% Save the dataXXX.mat files (with segmentations only) for safekeeping before processing
mkdir('segmentations_only')
%%  Process image objects
%
% BEFORE YOU RUN THIS, SAVE SEGMENTATION data00X.mat files TO ANOTHER FOLDER for safekeeping

%condList=[1:15];
condList=[16];

origDir=pwd;
sigma=0.5;
numLevels=3;

for i=1:length(condList)
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condList(i)})
    
    tools = improc2.launchImageObjectTools();
    
    % check if you have saved the segmentations-only .mat files
    if ~isfolder([Tcond.condName{rowTcond},'_segmentations']), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
    
    %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
    improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
    cd(parentDir)
end

%% exclude the bottom 2 slices
condList=1:45;
numBottomSlicesToExclude=2;

for i=1:length(condList)
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condList(i)})
    
    excludeSlices(numBottomSlicesToExclude)
    
    cd(parentDir)
end

%% Manually set threshold for smFISH (tmr) channel
% navigate to each condition folder before running this.
%improc2.launchThresholdGUI


%% save the processed (but not yet Gaussian-fitted) data to their own folder
%see save_dat_files.m
subfolderName='dataFiles2_processed';
if ~isdir(subfolderName)
    mkdir(subfolderName)
else
    error(sprintf("subfolder %s already exists, delete it first",subfolderName))
end
copyfile('*.mat',subfolderName)
%% set a uniform threshold for clampFISH channels
% this is done so that spot fitting knows how many spots to fit
%condList=2;% 2 started 117pm
%condList=1:45;
condList=[16:45,1:15];% started 124am 13Nov; about 12-13min to do thresholding. GFPser1 9-2 start 141am --> 20min until GFP_ser2 9-2. NoPrimser7 3-2 at 1124am. 10hr for 30+6. at 1232pm only on NoPrimser7 4-3 (7cells in 1hr). 2pm NoPrim_ser13 3-1. By 250pm it's done with fitting - doing extractions.
cd(parentDir)

%[200  200  200  200  200  200  200  200  200   200   200  200  200  200  200,... %NoPrimaries

%  1    2    3    4    5    6    7    8    9    10    11   12   13   14   15
thresholdList=... %This threshold list takes about 13hr to complete
[ 100  100  100  100  100  100  100  100  100   125   250  200  100   250  150,...
  450  150  400  250  250  300  200  200  400   500  1000  800  160  1000  650,... %GFP
 1900  270 1200  700  650 1000  850 1200 2000  3100  2900 1800  300  2300 1700,...%EGFR
 NaN  NaN  NaN]';


for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condID})
    
    threshold=thresholdList(condID);
    setThreshold('cya',threshold)
    
    cd(parentDir)
end

%%%%% %% Fit spots
cd(parentDir)

for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condID})
    
    threshold=thresholdList(condID);
    fitSpots('cya')
    
    cd(parentDir)
end


%%%%% %% Extract cell data --> Tcell

GFPratioLow=1.6;
GFPratioHigh=2;

TextraMeasurements=table('Size',[3,4],'VariableNames',{'columnName','area','valueToCompute','channel'},'VariableTypes',{'string','string','string','string'});
TextraMeasurements.columnName=          {  'avgNuclearGFP';    'avgCytoplasmicGFP'; 'avgCellGFP'};
TextraMeasurements.area=                {  'nuclear';          'cytoplasmic';       'cell'      };
TextraMeasurements.valueToCompute=      {  'average';          'average';           'average'   };
TextraMeasurements.channel=             {  'gfp';              'gfp';               'gfp'       };

cd(parentDir)

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    Tcell=extractCells(TextraMeasurements);
    Tcell.ratioGFPnuclearToCytoplasmic=Tcell.avgNuclearGFP./Tcell.avgCytoplasmicGFP;
    
    % asign predictedGFPstatus based on ratioGFPnuclearToCytoplasmic
    cellRows_with_ambiguous_GFP_status=find((all([Tcell.ratioGFPnuclearToCytoplasmic>GFPratioLow,Tcell.ratioGFPnuclearToCytoplasmic<GFPratioHigh],2)));
    if ~isempty(cellRows_with_ambiguous_GFP_status)
        warning('will remove cells in these rows because GFP status, since they are within the ambiguous range')
        display(cellRows_with_ambiguous_GFP_status)
    end
    Tcell.predictedGFPstatus=repmat({'ambiguous'},height(Tcell),1);
    Tcell.predictedGFPstatus(Tcell.ratioGFPnuclearToCytoplasmic<=GFPratioLow)={'negative'};
    Tcell.predictedGFPstatus(Tcell.ratioGFPnuclearToCytoplasmic>=GFPratioHigh)={'positive'};
    if any(strcmp(Tcell.predictedGFPstatus,'ambiguous'))
        warning('%i cells has ambiguous GFP status:',sum(strcmp(Tcell.predictedGFPstatus,'ambiguous')))
        Tcell(strcmp(Tcell.predictedGFPstatus,'ambiguous'),:)
    end

    % save Tcell in extracted directory
    cd(parentDir)
    extractedDataDirThisCond=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDirThisCond)
        mkdir(extractedDataDirThisCond)
    end
    writetable(Tcell,fullfile(extractedDataDirThisCond,'Tcell1.csv')); %for now, save it locally

end

%% Extract spot data --> Tspots

cd(parentDir)


numCellsWhosDataHasBeenExtracted=0;
tmrCutoffForExtraction=300;
cyaCutoffForExtraction=150;

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    % get spots
    Tspots=extractSpots({'tmr','cya'},'intensitiesAboveGivenCutoff',[tmrCutoffForExtraction,cyaCutoffForExtraction],'getFittedData',true);
    cd(parentDir)
    extractedDataDirThisCond=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDirThisCond)
        mkdir(extractedDataDirThisCond)
    end
    writetable(Tspots,fullfile(extractedDataDirThisCond,'Tspots1.csv')); %for now, save it locally
    
end
%% save the fitted data files
subfolderName='dataFiles3_fitted';
overwriteFiles=false;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condID})
    
    % removing extraneous directories
%     if isfolder('extractedData')
%         rmdir('extractedData','s')
%     end

    
    if overwriteFiles && isfolder(subfolderName)
        rmdir(subfolderName)
    end
    if ~isfolder(subfolderName)
        mkdir(subfolderName)
        copyfile('*.mat',subfolderName)
    else
        error(sprintf("delete directory %s or set overwriteFiles=true",subfolderName))
    end
    
    cd(parentDir)
end

% %% rename files (if needed)
% cd(parentDir)
% for i=1:length(condList)
%     condID=condList(i);
%     cd(Tcond.Extracted_60X_dir{Tcond.condID==condID})
%     
%     %Tcell=readtable('Tcell.csv');
%     %writetable(Tcell,'Tcell1.csv');
%     %delete('Tcell.csv')
%     %Tspots=readtable('Tspots.csv');
%     %writetable(Tspots,'Tspots1.csv');
%     %delete('Tspots.csv')
%     cd(parentDir)
%     
%     cd(parentDir)
% end
%% Remove ambiguous GFP status cells from Tcell and Tspots, add on cellID, then re-save --> Tspots2.csv, Tcell2.csv
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Extracted_60X_dir{Tcond.condID==condID})
    
    Tcell=readtable('Tcell1.csv');
    Tspots=readtable('Tspots1.csv');
    
    numAmbig=sum(strcmp(Tcell.predictedGFPstatus,'ambiguous'));
    
    Tcell=Tcell(~strcmp(Tcell.predictedGFPstatus,'ambiguous'),:); %remove ambiguous GFP status cells
    Tcell=[array2table([repmat(condID,height(Tcell),1),[1:height(Tcell)]'],'VariableNames',{'condID','cellID'}),Tcell];
    
    Tspots=innerjoin(Tspots,Tcell(:,{'cellID','arrayNum','objNum'}),'Keys',{'arrayNum','objNum'});
    Tspots=[array2table(repmat(condID,height(Tspots),1),'VariableNames',{'condID'}),Tspots(:,end:end),Tspots(:,1:end-1)]; %first two columns are condID, cellID
    
    writetable(Tcell,'Tcell2.csv');
    writetable(Tspots,'Tspots2.csv');
    
    cd(parentDir)
end

%% put all Tcell and Tspot table together into one
condList=1:45;

extractedDataDir='paper/extractedData/E143_screen/60X';
cd(parentDir)

TcellAll=table();
TspotsAll=table();

for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Extracted_60X_dir{Tcond.condID==condID})
    
    Tcell=readtable('Tcell2.csv');
    Tspots=readtable('Tspots2.csv');
    
    TcellAll=[TcellAll;Tcell];
    TspotsAll=[TspotsAll;Tspots];
    
    cd(parentDir)
end

writetable(TcellAll,fullfile(extractedDataDir,'TcellAll.csv'))
writetable(TspotsAll,fullfile(extractedDataDir,'TspotsAll.csv'))

%% Add onto Tcell the number of above-threshold conentional smFISH (tmr) spots
Tcell=readtable(fullfile(extractedDataDir,'TcellAll.csv'));
Tspots=readtable(fullfile(extractedDataDir,'TspotsAll.csv'));

tmrUseCellSpecificThreshold=true;
tmrThresholdUniform=1000;

% clampThresholdSweep=[1100]; % put multiple values to plot relationships % For EFFR_ser1, like cya threshold =4000 or 3500.

%if or(tmrThresholdUniform<tmrCutoffForExtraction,any(clampThresholdSweep<cyaCutoffForExtraction)), error("need to re-extract with lower cutoff if you want to use this threshold"),end
if tmrThresholdUniform<tmrCutoffForExtraction, error("need to re-extract tmr with lower cutoff if you want to use this tmr threshold"),end
% number of smFISH spots vs. number of clampFISH spots (with threshold sweep)
Tcell.tmrNumSpots=nan(height(Tcell),1);
%Tcell.cyaNumSpots=nan(height(Tcell),length(clampThresholdSweep));
Tcell.tmrCellSpecificThreshold=nan(height(Tcell),1);

condList=unique(Tcell.condID);
for iCond=1:length(condList)
    condID=condList(iCond)
    
    TcellThisCond=Tcell(Tcell.condID==condID,:);
    TspotsThisCond=Tspots(Tspots.condID==condID,:);
    cellIDsThisCond=sort(TcellThisCond.cellID,'ascend');
    assert(isequal(cellIDsThisCond,sort(unique(Tspots.cellID(Tspots.condID==condID)),'ascend')))
    
    for cellID=cellIDsThisCond'

        rowTcell=find(all([Tcell.condID==condID,Tcell.cellID==cellID],2));
        assert(length(rowTcell)==1)
        
        channel='tmr';
        tmrThresholdCellSpecific=unique(Tspots.threshold(all([Tspots.condID==condID,Tspots.cellID==cellID,strcmp(Tspots.channel,channel),Tspots.isGood],2)));
        if ~length(tmrThresholdCellSpecific)==1
            tmrThresholdCellSpecific=nan;
            warning('no tmr spots in Tspots for cellID=%i check the threshold you used to extract spots',cellID)
        end
        
        Tcell.tmrCellSpecificThreshold(rowTcell)=tmrThresholdCellSpecific;
        
        if tmrUseCellSpecificThreshold
            tmrThreshold=tmrThresholdCellSpecific;
        else
            tmrThreshold=tmrThresholdUniform;
        end
        tmrSpotRows=all([Tspots.condID==condID,Tspots.cellID==cellID,strcmp(Tspots.channel,channel),Tspots.isGood,Tspots.intensities>tmrThreshold],2);
        tmrNumSpots=sum(tmrSpotRows);
        Tcell.tmrNumSpots(rowTcell)=tmrNumSpots;
        
        %     channel='cya';
        %     for ii=1:length(clampThresholdSweep)
        %         threshold=clampThresholdSweep(ii);
        %         cyaSpotRows=all([Tspots.condID==condID,Tspots.cellID==cellID,strcmp(Tspots.channel,channel),Tspots.isGood,Tspots.intensities>threshold],2);
        %         cyaNumSpots=sum(cyaSpotRows);
        %         Tcell.cyaNumSpots(i,ii)=cyaNumSpots;
        %     end
    end
    
end

writetable(Tcell,fullfile(extractedDataDir,'TcellAll2_smFISHspots.csv'))
%% Add onto Tcell avg intensity of top N clampFISH spots, N=No. tmr spots
% not done for first 15 conditions, which are without primaries
Tcell=readtable(fullfile(extractedDataDir,'TcellAll2_smFISHspots.csv'));
Tspots=readtable(fullfile(extractedDataDir,'TspotsAll.csv'));

condList=16:45;

Tcell.intensityClampAtsmFISHcountAboveMin=nan(height(Tcell),1);
minsmFISHcount=20;
clampChannel='cya';
smFISHChannel='tmr';
for iCond=1:length(condList)
    condID=condList(iCond);
    TspotsThisCond=Tspots(Tspots.condID==condID,:);
    rowsOfTcell=find(Tcell.condID==condID)';
    for rowTcell=rowsOfTcell
        cellID=Tcell.cellID(rowTcell);
        smFISHnumSpots=Tcell.tmrNumSpots(rowTcell);
        
        if smFISHnumSpots>=minsmFISHcount
            clampSpotIntensities=sort(TspotsThisCond.intensities(all([TspotsThisCond.cellID==cellID,strcmp(TspotsThisCond.channel,clampChannel)],2)),'descend');
            avgIntensityTopNclampSpots=mean(clampSpotIntensities(1:smFISHnumSpots));
            Tcell.intensityClampAtsmFISHcountAboveMin(rowTcell)=avgIntensityTopNclampSpots;
        end
    end
end

writetable(Tcell,fullfile(extractedDataDir,'TcellAll3_smfishN_clampAvgIntNabove20.csv'))


%% check consistency between Tspots and Tcell
Tcell=readtable(fullfile(extractedDataDir,'TcellAll3_smfishN_clampAvgIntNabove20.csv'));
Tspots=readtable(fullfile(extractedDataDir,'TspotsAll.csv'));

% do they contain the same condIDs?
condIDList=sort(unique(Tcell.condID));
assert(isequal(condIDList,sort(unique(Tspots.condID))));

% are all rows in Tspots.isGood and Tcell.isGood true?
assert(all(Tspots.isGood))
assert(all(Tcell.isGood))

% do the condID and cellID all match?
assert(isequal(sortrows(Tcell(:,{'condID','cellID'})),sortrows(unique(Tspots(:,{'condID','cellID'}),'rows'))))

% are there any tmr spots below the user-selected threshold?
%assert(~any(Tspots.intensities(strcmp(Tspots.channel,'tmr'))<Tspots.threshold(strcmp(Tspots.channel,'tmr'))))

smFISHChannel='tmr';
%clampChannel='cya';

for condID=condIDList'
    
    TcellThisCond=Tcell(Tcell.condID==condID,:);
    TspotsThisCond=Tspots(Tspots.condID==condID,:);
    
    cellIDsThisCond=sort(TcellThisCond.cellID,'ascend');
    assert(isequal(cellIDsThisCond,sort(unique(Tspots.cellID(Tspots.condID==condID)),'ascend')))
    
    for cellID=cellIDsThisCond'
        
        % check tmr number of spots matches that in Tspots above threshold
        tmrNumSpotsFromTspots=sum(all([TspotsThisCond.cellID==cellID,strcmp(TspotsThisCond.channel,smFISHChannel),TspotsThisCond.intensities>=TspotsThisCond.threshold],2));
        tmrNumSpotsFromTcell=Tcell.tmrNumSpots(all([Tcell.condID==condID,Tcell.cellID==cellID],2));
        if~tmrNumSpotsFromTspots==tmrNumSpotsFromTcell
            error(sprintf('For condID=%i and cellID=%i, Tcell.tmrNumSpots is %i but Tspots has %i tmr spots\n',condID,cellID,tmrNumSpotsFromTcell,tmrNumSpotsFromTspots))
        end
        
    end
    
    
end


%% place top N clampFISH spot intensities for cells with N>=20 into TspotsClampTopN

condIDList=[16:45]';

TspotsClampTopN=table();

smFISHChannel='tmr';
clampChannel='cya';
minsmFISHcount=20;

for condID=condIDList'
    
    cellIDlistThisCond=Tcell.cellID(Tcell.condID==condID);
    
    for ii=1:length(cellIDlistThisCond)
        cellID=cellIDlistThisCond(ii);
        TspotsThisCell=Tspots(all([Tspots.condID==condID,Tspots.cellID==cellID],2),:);
        
        smFISHnumSpots=sum(all([strcmp(TspotsThisCell.channel,smFISHChannel),TspotsThisCell.intensities>=TspotsThisCell.threshold],2));
        %assert(smFISHnumSpots==Tcell.tmrNumSpots(all([Tcell.condID==condID,Tcell.cellID==cellID],2)));
        %%Added 12Feb2021
        
        if smFISHnumSpots>=minsmFISHcount
            clampSpotIntensities=sort(TspotsThisCell.intensities(strcmp(TspotsThisCell.channel,clampChannel)),'descend');
            lowestClampIntensityInTopN=clampSpotIntensities(smFISHnumSpots);
            
            TspotsClampTopN_new=TspotsThisCell(all([strcmp(TspotsThisCell.channel,clampChannel),TspotsThisCell.intensities>=lowestClampIntensityInTopN],2),:);
            
            if height(TspotsClampTopN)==0
                TspotsClampTopN=TspotsClampTopN_new;
            else
                TspotsClampTopN=[TspotsClampTopN;TspotsClampTopN_new];
            end
            TspotsClampTopN_new=[];
            
            
        end
    end
    
    
end
writetable(TspotsClampTopN,fullfile(extractedDataDir,filesep,'TspotsClampTopN.csv'))


%% Tcond

%
% for i=1:length(condUniqueList)
%     condID=condUniqueList(i);
%
%     avgIntensityTopNclamp=mean(TspotsClampTopN.intensities(all([TspotsClampTopN.condID==condID,strcmp(TspotsClampTopN.channel,'cya')],2)));
%     Tcond.avgIntensityTopNclamp(Tcond.condID==condID)=avgIntensityTopNclamp;
%
%     stdIntensityTopNclamp=std(TspotsClampTopN.intensities(all([TspotsClampTopN.condID==condID,strcmp(TspotsClampTopN.channel,'cya')],2)));
%     Tcond.stdIntensityTopNclamp(Tcond.condID==condID)=stdIntensityTopNclamp;
%
%     numSpotTopN=sum(all([TspotsClampTopN.condID==condID,strcmp(TspotsClampTopN.channel,'cya')],2));
%     Tcond.numSpotTopN(Tcond.condID==condID)=numSpotTopN;
% end
%
%

%% Add summary metrics to TcondResults
condIDList=[16:45]';

TspotsClampTopN=readtable(fullfile(extractedDataDir,filesep,'TspotsClampTopN.csv'));

TcondResults=Tcond(:,1:4);
TcondResults=[TcondResults,array2table(nan(height(TcondResults),10),'VariableNames',{'Q0','Q10','Q25','median','Q75','Q90','mean','std','Nspots','Ncells'}')];

% find on-target (topN clampFISH) statistics:
for condID=condIDList'
    
    rowTcondResults=find(TcondResults.condID==condID); assert(length(rowTcondResults)==1)
    spotList=TspotsClampTopN.intensities(TspotsClampTopN.condID==condID);
    
    TcondResults.Q0(rowTcondResults)=prctile(spotList,0);
    TcondResults.Q10(rowTcondResults)=prctile(spotList,10);
    TcondResults.Q25(rowTcondResults)=prctile(spotList,25);
    TcondResults.median(rowTcondResults)=prctile(spotList,50);
    TcondResults.Q75(rowTcondResults)=prctile(spotList,75);
    TcondResults.Q90(rowTcondResults)=prctile(spotList,90);
    
    TcondResults.mean(rowTcondResults)=mean(spotList);
    TcondResults.std(rowTcondResults)=std(spotList);
    
    TcondResults.Nspots(rowTcondResults)=length(spotList);
    TcondResults.Ncells(rowTcondResults)=length(unique(TspotsClampTopN.cellID(TspotsClampTopN.condID==condID)));
    
end

TcondResults.pctStd=100*TcondResults.std./TcondResults.mean;
TcondResults

writetable(TcondResults,fullfile(extractedDataDir,filesep,'TcondResults.csv'))


%% write this condition's Tcell & Tspots to csv file
% writetable(Tcell,'Tcell.csv')
% writetable(Tspots,'Tspots.csv')
% 
% 
% %%
% Tcell=readtable('Tcell.csv');
% Tspots=readtable('Tspots.csv');

