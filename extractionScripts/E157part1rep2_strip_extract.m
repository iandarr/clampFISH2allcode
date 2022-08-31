% E157part1rep2_strip_extract
% readout stripping
%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E157part1rep2_strip_Tcond.csv'));
Tscan=readtable(fullfile(parentDir,filesep,'paper/experiments/E157part1rep2_strip_Tscan.xlsx'));
subregionDir='SubregionArray2x2/Subregion_1_r1_c1';
extractedDataDirForAllCond='paper/extractedData/E157part1rep2_strip/';
scanIDlist=Tscan.scanID';
%% Pre-extraction steps
% for each row in Tscan, only need to use first subregion's stitch (Subregion1_r1_c2)
% after stitching, append a 001 to the end of the filenames
% copy-paste R1_DAPI_50ms.tif as dapi.tif

%% Segmentation
%navigate to a subfolder and run this. Segment cells.

%improc2.segmentGUI.SegmentGUI()

%% save the segmentations to a subfolder for safekeeping
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
warning('short scanIDlist')
cd(parentDir)
for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    rowTscan=find(Tscan.scanID==scanID);
    cd([Tscan.scanDir{rowTscan},filesep,subregionDir])
 
    if overwriteFiles && isfolder(subfolderName)
        rmdir(subfolderName,'s')
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
display('processing image objects')
origDir=pwd;
sigma=0.4;
numLevels=3;

cd(parentDir)
for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    rowTscan=find(Tscan.scanID==scanID);
    cd([Tscan.scanDir{rowTscan},filesep,subregionDir])
    
    tools = improc2.launchImageObjectTools();
    
    % check if you have saved the segmentations-only .mat files
    if ~isfolder('dataFiles1_segmentations'), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
    
    %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
    improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
    cd(parentDir)
end

%% save the processed data files for safekeeping
subfolderName='dataFiles2_processed';
overwriteFiles=true;

cd(parentDir)
for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    rowTscan=find(Tscan.scanID==scanID);
    cd([Tscan.scanDir{rowTscan},filesep,subregionDir])
    
    if overwriteFiles && isfolder(subfolderName)
        rmdir(subfolderName,'s')
    end
    if ~isfolder(subfolderName)
        mkdir(subfolderName)
        copyfile('*.mat',subfolderName)
    else
        error(sprintf("delete directory %s or set overwriteFiles=true",subfolderName))
    end

    cd(parentDir)
end

%% Extract cell data --> Tcell

display('extracting cell data')
cd(parentDir)

TextraMeasurements=table('Size',[4,4],'VariableNames',{'columnName','area','valueToCompute','channel'},'VariableTypes',{'string','string','string','string'});
TextraMeasurements.columnName=          {  'avgCellYFPc';      'avgNoncellYFPc'; 'avgNuclearDAPI';  'avgCytoplasmDAPI'};
TextraMeasurements.area=                {  'cell';             'nonCell';        'nuclear'       ;  'cytoplasmic'      };
TextraMeasurements.valueToCompute=      {  'average';          'average';        'average'       ;  'average'         };
TextraMeasurements.channel=             {  'gfpc';             'gfpc';           'dapi'          ;  'dapi'            };

cd(parentDir)
for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    rowTscan=find(Tscan.scanID==scanID);
    cd([Tscan.scanDir{rowTscan},filesep,subregionDir])
    
    Tcell=extractCells();
    %Tcell=extractCells(TextraMeasurements);
    
    A488_gene=Tscan.A488_gene(rowTscan);
    A555_gene=Tscan.A555_gene(rowTscan);
    A594_gene=Tscan.A594_gene(rowTscan);
    A647_gene=Tscan.A647_gene(rowTscan);
    
    nCell=height(Tcell);
    Tcell=[Tcell,array2table([repmat(A488_gene,nCell,1),repmat(A555_gene,nCell,1),repmat(A594_gene,nCell,1),repmat(A647_gene,nCell,1)],'VariableNames',{'A488_gene','A555_gene','A594_gene','A647_gene'})];
    % save Tcell in extracted directory
    cd(parentDir)
    extractedDataDir=Tscan.extractedDataDir{rowTscan};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    
    writetable(Tcell,fullfile(extractedDataDir,'Tcell1.csv')); %for now, save it locally
   
end

%% Extract spot data --> Tspots
% NOTE: need to add the gaussian fitting before this, then also extract those
cd(parentDir)
%scanIDlist=1;


channelsToExtract={'all'};
intensityCutoffs=30;%repmat(30,1,length(channelsToExtract));

for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    rowTscan=find(Tscan.scanID==scanID);
    cd([Tscan.scanDir{rowTscan},filesep,subregionDir])
    
    % get spots
    Tspots=extractSpots(channelsToExtract,'intensitiesAboveGivenCutoff',intensityCutoffs,'getFittedData',false);
    cd(parentDir)
    extractedDataDir=Tscan.extractedDataDir{rowTscan};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    
end

%% put all Tcell and Tspot table together into one
%scanIDlist=1;

if ~isfolder(extractedDataDirForAllCond)
        mkdir(extractedDataDirForAllCond)
end
cd(parentDir)

TcellAll=table();
TspotsAll=table();

for i=1:length(scanIDlist)
    scanID=scanIDlist(i);
    
    rowTscan=find(Tscan.scanID==scanID);
    cellType=Tscan.cellType{rowTscan};
    extractedDataDir=Tscan.extractedDataDir{rowTscan};
    cd(extractedDataDir)
    
    Tcell=readtable('Tcell1.csv');
    Tspots=readtable('Tspots1.csv');
    
    TcellAll=[TcellAll;[array2table(repmat(scanID,height(Tcell),1),'VariableNames',{'scanID'}),array2table(repmat({cellType},height(Tcell),1),'VariableNames',{'cellType'}),Tcell]];
    TspotsAll=[TspotsAll;[array2table(repmat(scanID,height(Tspots),1),'VariableNames',{'scanID'}),Tspots]];
    
    cd(parentDir)
end

writetable(TcellAll,fullfile(extractedDataDirForAllCond,'Tcell1All.csv'))
writetable(TspotsAll,fullfile(extractedDataDirForAllCond,'Tspots1All.csv'))

%% import data

TspotsAll=readtable(fullfile(extractedDataDirForAllCond,'Tspots1All.csv'));
TcellAll1=readtable(fullfile(extractedDataDirForAllCond,'Tcell1All.csv'));

%%
Tthresh=readtable(fullfile(extractedDataDirForAllCond,'Tthresholds.xlsx'));
temp1=cell(height(Tthresh),1);
for i=1:height(Tthresh)
   temp1{i}=str2num(Tthresh.scanIDs{i});    
end
Tthresh.scanIDs=temp1;
%% make TcellAll2,which has numSpots data (spot intensities above threshold) for before & after stripping
% first 
TcellPlotPaired=table();
condID_readList=Tcond.condID(all([ismember(Tcond.condID,1:10),startsWith(Tcond.readoutCycle,'cycle')],2))';
%
TcellAll2=TcellAll1(:,[2:4,1,5:end]);
nCell=height(TcellAll2);
TcellAll2=[TcellAll2(:,1:end-3),...
                     array2table(cell(nCell,2),'VariableNames',{'read_A488_channel','strip_A488_channel'}),array2table(nan(nCell,1),'VariableNames',{'A488_thresh'}),array2table(nan(nCell,2),'VariableNames',{'read_A488_numSpots','strip_A488_numSpots'}),...
           TcellAll2(:,end-2),...
                     array2table(cell(nCell,2),'VariableNames',{'read_A555_channel','strip_A555_channel'}),array2table(nan(nCell,1),'VariableNames',{'A555_thresh'}),array2table(nan(nCell,2),'VariableNames',{'read_A555_numSpots','strip_A555_numSpots'}),...
           TcellAll2(:,end-1),...
                     array2table(cell(nCell,2),'VariableNames',{'read_A594_channel','strip_A594_channel'}),array2table(nan(nCell,1),'VariableNames',{'A594_thresh'}),array2table(nan(nCell,2),'VariableNames',{'read_A594_numSpots','strip_A594_numSpots'}),...
           TcellAll2(:,end),...
                     array2table(cell(nCell,2),'VariableNames',{'read_A647_channel','strip_A647_channel'}),array2table(nan(nCell,1),'VariableNames',{'A647_thresh'}),array2table(nan(nCell,2),'VariableNames',{'read_A647_numSpots','strip_A647_numSpots'})];
for i=1:height(TcellAll2)
    scanID=TcellAll2.scanID(i);
    arrayNum=TcellAll2.arrayNum(i);
    objNum=TcellAll2.objNum(i);
    
    prefixList={'A488','A555','A594','A647'};
    for ii=1:length(prefixList)
        prefix=prefixList{ii};
        gene=TcellAll2.([prefix,'_gene']){i};
        % A488 gene
        indTthresh=strcmp(Tthresh.gene,gene);
        assert(sum(indTthresh)==1)
        
        if ismember(scanID,Tthresh.scanIDs{indTthresh})
            read_channelRajlabimagetools=Tthresh.read_channelRajlabimagetools(indTthresh);
            strip_channelRajlabimagetools=Tthresh.strip_channelRajlabimagetools(indTthresh);
            threshold=Tthresh.threshold(indTthresh);
            TcellAll2.(['read_',prefix,'_channel'])(i)=read_channelRajlabimagetools;
            TcellAll2.(['strip_',prefix,'_channel'])(i)=strip_channelRajlabimagetools;
            TcellAll2.([prefix,'_thresh'])(i)=threshold;
            
            % get numSpots
            TspotsThisCellThisChannel_read=TspotsAll(all([TspotsAll.scanID==scanID,TspotsAll.arrayNum==arrayNum,TspotsAll.objNum==objNum,strcmp(TspotsAll.channel,read_channelRajlabimagetools),TspotsAll.intensities>=threshold],2),:);
            numSpots_read=height(TspotsThisCellThisChannel_read);
            TcellAll2.(['read_',prefix,'_numSpots'])(i)=numSpots_read;

            TspotsThisCellThisChannel_strip=TspotsAll(all([TspotsAll.scanID==scanID,TspotsAll.arrayNum==arrayNum,TspotsAll.objNum==objNum,strcmp(TspotsAll.channel,strip_channelRajlabimagetools),TspotsAll.intensities>=threshold],2),:);
            numSpots_strip=height(TspotsThisCellThisChannel_strip);
            TcellAll2.(['strip_',prefix,'_numSpots'])(i)=numSpots_strip;

        end
    end
    
    
end
%
writetable(TcellAll2,fullfile(extractedDataDirForAllCond,'TcellAll2.csv'))