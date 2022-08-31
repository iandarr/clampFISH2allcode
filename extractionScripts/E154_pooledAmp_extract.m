% E154_pooledAmp_extract
% pooled amplification

%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E154_pooled_amp_Tcond.xlsx');


%% Segmentation
%navigate to a subfolder and run this. Segment cells.

improc2.segmentGUI.SegmentGUI

%% save the segmentations to a subfolder for safekeeping
subfolderName='dataFiles1_segmentations';
overwriteFiles=false;
condList=25:44;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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

%%  Process image objects
% BEFORE YOU RUN THIS, SAVE SEGMENTATIONS TO ANOTHER FOLDER
condList=25:44;
origDir=pwd;
sigma=0.5;
numLevels=3;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    tools = improc2.launchImageObjectTools();
    
    % check if you have saved the segmentations-only .mat files
    if ~isfolder('dataFiles1_segmentations'), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
    
    %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
    improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
    cd(parentDir)
end
%% save the processed data files for safekeeping
subfolderName='dataFiles2_processed';
overwriteFiles=false;
condList=25:44;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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

%% Thresholding (did not do)
improc2.launchThresholdGUI
%% Extract cell data --> Tcell
% left out extra measurements since not all have nuclear mask
cd(parentDir)

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    Tcell=extractCells();
    
    % save Tcell in extracted directory
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tcell,fullfile(extractedDataDir,'Tcell1.csv')); %for now, save it locally

end

%% Extract spot data --> Tspots
% NOTE: need to add the gaussian fitting before this, then also extract those
cd(parentDir)

tmrCutoffForExtraction=300;
cyaCutoffForExtraction=150;

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    % get spots
    Tspots=extractSpots({'tmr','cya'},'intensitiesAboveGivenCutoff',[tmrCutoffForExtraction,cyaCutoffForExtraction],'getFittedData',false);
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    
end

%% To Tspots1.csv and Tcell1.csv: Add on cellID, then re-save --> Tspots2.csv, Tcell2.csv
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Extracted_60X_dir{rowTcond})
    
    Tcell=readtable('Tcell1.csv');
    Tspots=readtable('Tspots1.csv');
        
    Tcell=[array2table([repmat(condID,height(Tcell),1),[1:height(Tcell)]'],'VariableNames',{'condID','cellID'}),Tcell];
    
    Tspots=innerjoin(Tspots,Tcell(:,{'cellID','arrayNum','objNum'}),'Keys',{'arrayNum','objNum'});
    Tspots=[array2table(repmat(condID,height(Tspots),1),'VariableNames',{'condID'}),Tspots(:,end:end),Tspots(:,1:end-1)]; %first two columns are condID, cellID
    
    writetable(Tcell,'Tcell2.csv');
    writetable(Tspots,'Tspots2.csv');
    
    cd(parentDir)
end


%% sub-sample 40 cells per condition --> save in Tspots3.csv, Tcell3.csv
cd(parentDir)
extractedDataDir='paper/extractedData/E154_pooled_amp/s6 60X p2';

condList=25:44;
numCellsToSubsample=40; % cells to sub-sample

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Extracted_60X_dir{rowTcond})
    
    Tcell=readtable('Tcell2.csv');
    
    cellIDlist=unique(Tcell.cellID);
    if length(cellIDlist)<numCellsToSubsample
        error('attempting to subsample %i cells but only have %i with condID=%i',numCellsToSubsample,length(cellIDlist),condID)
    end
    
    rng(0)
    cellIDselected=sort(cellIDlist(randperm(length(cellIDlist),numCellsToSubsample)));
    
    % subsample Tcell
   indSelectedTcell=ismember(Tcell.cellID,cellIDselected);
   writetable(Tcell(indSelectedTcell,:),'Tcell3.csv');
   
   % Tspots_subsample
   Tspots=readtable('Tspots2.csv');
   indSelectedTspots=ismember(Tspots.cellID,cellIDselected);
   writetable(Tspots(indSelectedTspots,:),'Tspots3_1.csv');
   cd(parentDir)
end

%% put all Tcell and Tspot table together into one

extractedDataDir='paper/extractedData/E154_pooled_amp/s6 60X p2';
cd(parentDir)

TcellAll=table();
TspotsAll=table();

for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Extracted_60X_dir{Tcond.condID==condID})
    
    Tcell=readtable('Tcell3.csv');
    Tspots=readtable('Tspots3_1.csv');
    
    TcellAll=[TcellAll;Tcell];
    TspotsAll=[TspotsAll;Tspots];
    
    cd(parentDir)
end

writetable(TcellAll,fullfile(extractedDataDir,'TcellAll_3.csv'))
writetable(TspotsAll,fullfile(extractedDataDir,'TspotsAll_3_1.csv'))
%% put top K brightest clampFISH spots into a table TspotsClampTopK
cd(parentDir)
Tspots=readtable(fullfile(extractedDataDir,'TspotsAll_3_1.csv'));

K=10000;
TspotsClampTopN=table();
condList=25:44;
clampChannel='cya';
TspotsClampTopK=table();
for condID=condList
    
	TspotsClampThisCond=Tspots(all([Tspots.condID==condID,strcmp(Tspots.channel,clampChannel)],2),:);
    TspotsClampThisCond=sortrows(TspotsClampThisCond,{'intensities'},'descend');
    %condID
    %numSpots=height(TspotsClampThisCond)
    TspotsClampThisCond=TspotsClampThisCond(1:K,:);
    
    TspotsClampTopK=[TspotsClampTopK;TspotsClampThisCond(1:K,:)];
    
end
    
writetable(TspotsClampTopK,fullfile(extractedDataDir,'TspotsClampTopK3_1.csv'))

%% Add summary metrics to TcondResults
TspotsClampTopK=readtable(fullfile(extractedDataDir,'TspotsClampTopK3_1.csv'));

condList=25:44;

TcondResults=Tcond(:,1:10);
TcondResults=[TcondResults,array2table(nan(height(TcondResults),8),'VariableNames',{'Q0','Q10','Q25','median','Q75','Q90','mean','std'}')];

% find on-target (topN clampFISH) statistics:
for condID=condList
    
    rowTcondResults=find(TcondResults.condID==condID); assert(length(rowTcondResults)==1)
    spotList=TspotsClampTopK.intensities(TspotsClampTopK.condID==condID);
    
    TcondResults.Q0(rowTcondResults)=prctile(spotList,0);
    TcondResults.Q10(rowTcondResults)=prctile(spotList,10);
    TcondResults.Q25(rowTcondResults)=prctile(spotList,25);
    TcondResults.median(rowTcondResults)=prctile(spotList,50);
    TcondResults.Q75(rowTcondResults)=prctile(spotList,75);
    TcondResults.Q90(rowTcondResults)=prctile(spotList,90);
    
    TcondResults.mean(rowTcondResults)=mean(spotList);
    TcondResults.std(rowTcondResults)=std(spotList);
    
end

TcondResults.pctStd=100*TcondResults.std./TcondResults.mean;


% Add on % difference between means, medians, Q10
condList_Pool=25:34;
condList_Alone=35:44;
%condList_Pool=TcondResults.condID(all([ismember(TcondResults.condID,condList),strcmp(TcondResults.amplifiersIncluded,'pool')],2))';

TcondResults.MeanPctDiff=nan(height(TcondResults),1);
TcondResults.MedianPctDiff=nan(height(TcondResults),1);
TcondResults.Q10PctDiff=nan(height(TcondResults),1);
for condID=condList_Pool
    
    
    rowTcondResults_pool=find(TcondResults.condID==condID);
    series=TcondResults.amplifierSeries(rowTcondResults_pool);
    rowTcondResults_alone=find(all([TcondResults.amplifierSeries==series,strcmp(TcondResults.amplifiersIncluded,'alone'),ismember(TcondResults.condID,condList_Alone)],2));
    
    TcondResults.MeanPctDiff(rowTcondResults_pool)=100*(TcondResults.mean(rowTcondResults_pool)/TcondResults.mean(rowTcondResults_alone) - 1);
    TcondResults.MedianPctDiff(rowTcondResults_pool)=100*(TcondResults.median(rowTcondResults_pool)/TcondResults.median(rowTcondResults_alone) - 1);
    TcondResults.Q10PctDiff(rowTcondResults_pool)=100*(TcondResults.Q10(rowTcondResults_pool)/TcondResults.Q10(rowTcondResults_alone) - 1);

end
TcondResults
writetable(TcondResults,fullfile(extractedDataDir,'TcondResults_1.csv'))

%% set a uniform threshold (single threshold per condition) for clampFISH channels
% this is done so that spot fitting knows how many regional maxima to fit. If you fit
% all the very dim regional maxima, it will take a very long time.

condList=25:44;
cd(parentDir)

% thresholds for condID=
%    25    26    27    28    29    30    31    32    33    34    35    36    37    38    39    40    41    42    43    44
thresholdList=...
[1900;	900;	400;	1100;	500;	1400;	1500;	2400;	2400;	1500;	2300;	800;	400;	900;	900;	2200;	2300;	1800;	1600;	1600];
% these numbers are round(0.9 * intensity of 10,000'th brightest regional maximum of the condition)

for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condID})
    
    threshold=thresholdList(i);
    setThreshold('cya',threshold)
    
    cd(parentDir)
end

%% save the thresholded data files for safekeeping
subfolderName='dataFiles3_thresholded';
overwriteFiles=false;
condList=25:44;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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


%% Fit spots %conID 25 at 23:19pm; condID=28 start 11:34pm (5min per condition) finish at 1am; condID=35 starting at 12:07am
cd(parentDir)

for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_60X_dir{Tcond.condID==condID})
    fprintf('condID=%i at %s',condID,char(datetime))
    
    fitSpots('cya')
    
    cd(parentDir)
end
%% save the fitted data files for safekeeping
subfolderName='dataFiles4_fitted';
overwriteFiles=false;
condList=25:44;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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

%% Extract spot data, also with fitted intensities --> Tspots3_2
cd(parentDir)

tmrCutoffForExtraction=300;
cyaCutoffForExtraction=150;

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    % get spots
    Tspots=extractSpots({'tmr','cya'},'intensitiesAboveGivenCutoff',[tmrCutoffForExtraction,cyaCutoffForExtraction],'getFittedData',true);
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tspots,fullfile(extractedDataDir,'Tspots3_2.csv')); %for now, save it locally
    
end

%% To Tspots3_2.csv, add on cellID, then re-save --> Tspots3_3.csv
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Extracted_60X_dir{rowTcond})
    
    Tcell=readtable('Tcell3.csv');% hass cellID
    Tspots=readtable('Tspots3_2.csv');
    
    Tspots=innerjoin(Tspots,Tcell(:,{'cellID','arrayNum','objNum'}),'Keys',{'arrayNum','objNum'});
    Tspots=[array2table(repmat(condID,height(Tspots),1),'VariableNames',{'condID'}),Tspots(:,end:end),Tspots(:,1:end-1)]; %first two columns are condID, cellID
    
    writetable(Tspots,'Tspots3_3.csv');
    
    cd(parentDir)
end
%% put Tspot3_3.csv tables (with fitted spot intensities) together into one

extractedDataDir='extractedData/E154_pooled_amp/s6 60X p2';
cd(parentDir)
TspotsAll=table();
for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Extracted_60X_dir{Tcond.condID==condID})
    
    Tspots=readtable('Tspots3_3.csv');
    
    TspotsAll=[TspotsAll;Tspots];
    
    cd(parentDir)
end
writetable(TspotsAll,fullfile(extractedDataDir,'TspotsAll_3_3.csv'))
%% put top K brightest (by unfitted intensity) clampFISH spots into a table TspotsClampTopK
% now the data have a fitted spot intensity associated with it as well
cd(parentDir)
Tspots=readtable(fullfile(extractedDataDir,'TspotsAll_3_3.csv'));

K=10000;
TspotsClampTopN=table();
condList=25:44;
clampChannel='cya';
TspotsClampTopK=table();
for condID=condList
    
	TspotsClampThisCond=Tspots(all([Tspots.condID==condID,strcmp(Tspots.channel,clampChannel)],2),:);
    TspotsClampThisCond=sortrows(TspotsClampThisCond,{'intensities'},'descend');
    %condID
    %numSpots=height(TspotsClampThisCond)
    TspotsClampThisCond=TspotsClampThisCond(1:K,:);
    
    TspotsClampTopK=[TspotsClampTopK;TspotsClampThisCond(1:K,:)];
    
end
    
writetable(TspotsClampTopK,fullfile(extractedDataDir,'TspotsClampTopK3_3.csv'))


%% Add summary metrics to TcondResults
TspotsClampTopK=readtable(fullfile(extractedDataDir,'TspotsClampTopK3_3.csv'));

condList=25:44;

TcondResults=Tcond(:,1:10);
TcondResults=[TcondResults,array2table(nan(height(TcondResults),8),'VariableNames',{'Q0','Q10','Q25','median','Q75','Q90','mean','std'}')];

% find on-target (topN clampFISH) statistics:
for condID=condList
    
    rowTcondResults=find(TcondResults.condID==condID); assert(length(rowTcondResults)==1)
    spotList=TspotsClampTopK.amplitudeFitted(TspotsClampTopK.condID==condID);
    
    TcondResults.Q0(rowTcondResults)=prctile(spotList,0);
    TcondResults.Q10(rowTcondResults)=prctile(spotList,10);
    TcondResults.Q25(rowTcondResults)=prctile(spotList,25);
    TcondResults.median(rowTcondResults)=prctile(spotList,50);
    TcondResults.Q75(rowTcondResults)=prctile(spotList,75);
    TcondResults.Q90(rowTcondResults)=prctile(spotList,90);
    
    TcondResults.mean(rowTcondResults)=mean(spotList);
    TcondResults.std(rowTcondResults)=std(spotList);
    
end

TcondResults.pctStd=100*TcondResults.std./TcondResults.mean;
TcondResults

% Add on % difference between means, medians, Q10
condList_Pool=25:34;
condList_Alone=35:44;
%condList_Pool=TcondResults.condID(all([ismember(TcondResults.condID,condList),strcmp(TcondResults.amplifiersIncluded,'pool')],2))';

TcondResults.MeanPctDiff=nan(height(TcondResults),1);
TcondResults.MedianPctDiff=nan(height(TcondResults),1);
TcondResults.Q10PctDiff=nan(height(TcondResults),1);
for condID=condList_Pool
    
    
    rowTcondResults_pool=find(TcondResults.condID==condID);
    series=TcondResults.amplifierSeries(rowTcondResults_pool);
    rowTcondResults_alone=find(all([TcondResults.amplifierSeries==series,strcmp(TcondResults.amplifiersIncluded,'alone'),ismember(TcondResults.condID,condList_Alone)],2));
    
    TcondResults.MeanPctDiff(rowTcondResults_pool)=100*(TcondResults.mean(rowTcondResults_pool)/TcondResults.mean(rowTcondResults_alone) - 1);
    TcondResults.MedianPctDiff(rowTcondResults_pool)=100*(TcondResults.median(rowTcondResults_pool)/TcondResults.median(rowTcondResults_alone) - 1);
    TcondResults.Q10PctDiff(rowTcondResults_pool)=100*(TcondResults.Q10(rowTcondResults_pool)/TcondResults.Q10(rowTcondResults_alone) - 1);

end
TcondResults
writetable(TcondResults,fullfile(extractedDataDir,'TcondResults_Fitted.csv'))


% %% Plot the relationship between cell area and spots
% Tcell_subsample.numClampSpotsTopK=nan(height(Tcell_subsample),1);
% 
% for i=1:height(Tcell_subsample)
%     condID=Tcell_subsample.condID(i);
%     cellID=Tcell_subsample.cellID(i);
%     numClampSpotsTopK=sum(all([TspotsClampTopK.condID==condID,TspotsClampTopK.cellID==cellID],2));
%     Tcell_subsample.numClampSpotsTopK(i)=numClampSpotsTopK;
% end
% 
% X=Tcell_subsample.areaCell;
% Y=Tcell_subsample.numClampSpotsTopK;
% G=cellstr(num2str(Tcell_subsample.condID));
% gscatter(X,Y,G)
%% boxplot of cellArea
% boxplot(X,G)
% 
% GroupOrderNum=[35 25 36 26 37 27 38 28 39 29 40 30 41 31 42 32 43 33 44 34];
% GroupOrder=strsplit(num2str(GroupOrderNum));
% boxplot(X,G,'GroupOrder',GroupOrder)
% ax=gca;
% %
% Tdum1=array2table(GroupOrderNum','VariableNames',{'condID'})
% Tdum2=Tcond(ismember(Tcond.condID,GroupOrderNum),{'condID','amplifierSeries','amplifiersIncluded'});
% Tdum2.amplifierSeries=replace(cellstr(num2str(Tdum2.amplifierSeries)),' ','');
% Tdum2.amplifiersIncluded=replace(Tdum2.amplifiersIncluded,{'alone','pool'},{'a','p'});
% %
% Tdum2.label=join([Tdum2.amplifierSeries,Tdum2.amplifiersIncluded],'');
% Tdum3=join(Tdum1,Tdum2)
% XTickLabel=Tdum3.label;
% ax.XTickLabel=XTickLabel;
