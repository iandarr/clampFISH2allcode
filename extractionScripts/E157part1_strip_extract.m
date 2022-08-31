% E157part1_strip_extract
% readout stripping
%% Table of all conditions

%parentDir='/Volumes/IAND_04';
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E157_validation_Tcond.xlsx')

%% Segmentation
%navigate to a subfolder and run this. Segment cells.

%improc2.segmentGUI.SegmentGUI(11)

%% copy data0XX.mat files (with segmentations) to post-strip folders
cd(parentDir)

moveCondIDs=[2,4,6,8,10];
for moveCondID=moveCondIDs %:length(fixedCondIDs)
    moveCondID
    fixedCondID=Tcond.ReferenceCondID(moveCondID)
    fixedTiffDir=Tcond.Tiff_20X_dir{fixedCondID};
    moveTiffDir=Tcond.Tiff_20X_dir{moveCondID};
        
    dataFileNames=listFilesWithExtension('.mat','directoryPath',fixedTiffDir);
    dataFileNames=dataFileNames(~startsWith(dataFileNames,'.'));
    
    for ii=1:length(dataFileNames)
        
        %%%% Per data0XX.mat file %%%%
        % copy-paste data file 
        dataFileName=dataFileNames{ii};
        fixedDataFilePath=fullfile(fixedTiffDir,dataFileName);
        copyfile(fixedDataFilePath,moveTiffDir);

        % change raw image file reference
        cd(moveTiffDir)
        changeFilePathInMatFile(dataFileName,moveTiffDir)
        cd(parentDir)

    end
    %%% For all objects in the folder, apply an XY shift to each non-reference object's bounding box
    moveOffsetX=Tcond.Offset20X_x(moveCondID);
    moveOffsetY=Tcond.Offset20X_y(moveCondID);
    
    cd(moveTiffDir)
    tools=improc2.launchImageObjectTools;
    while tools.iterator.continueIteration
        objectHandle=tools.objectHandle;
        shiftCurrentObjSegmentationMask(objectHandle,moveOffsetX,moveOffsetY); % main function that shifts the segmentations
        tools.iterator.goToNextObject
    end
        
    cd(parentDir)
end

cd(moveTiffDir)
%% Check the XY offsets by opening segmentation GUI on both 'fixed' and 'move' datasets
% cd(parentDir)
% moveCondIDs=[4];
% for moveCondID=moveCondIDs %:length(fixedCondIDs)
%     moveCondID
%     fixedCondID=Tcond.ReferenceCondID(moveCondID)
%     fixedTiffDir=Tcond.Tiff_20X_dir{fixedCondID};
%     moveTiffDir=Tcond.Tiff_20X_dir{moveCondID};
%     % open segmentGUI for fixed images
%     cd(fixedTiffDir)
%     improc2.segmentGUI.SegmentGUI
%     cd(parentDir)
%     % open segmentGUI for move images
%     cd(moveTiffDir)
%     improc2.segmentGUI.SegmentGUI
%     cd(parentDir)
%     
% end

%% save the segmentations to a subfolder for safekeeping
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
condList=1:10;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_20X_dir{rowTcond})
 
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
condList=1:10;
display('processing image objects')
origDir=pwd;
sigma=0.4;
numLevels=3;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_20X_dir{rowTcond})
    
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
condList=1:10;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_20X_dir{rowTcond})
 
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

% %% Thresholding (did not do)
%improc2.launchThresholdGUI

%% Extract cell data --> Tcell
% left out extra measurements since not all have nuclear mask
condList=1:10;
display('extracting cell data')
cd(parentDir)

TextraMeasurements=table('Size',[4,4],'VariableNames',{'columnName','area','valueToCompute','channel'},'VariableTypes',{'string','string','string','string'});
TextraMeasurements.columnName=          {  'avgCellYFPc';      'avgNoncellYFPc'; 'avgNuclearDAPI';  'avgCytoplasmDAPI'};
TextraMeasurements.area=                {  'cell';             'nonCell';        'nuclear'       ;  'cytoplasmic'      };
TextraMeasurements.valueToCompute=      {  'average';          'average';        'average'       ;  'average'         };
TextraMeasurements.channel=             {  'gfpc';             'gfpc';           'dapi'          ;  'dapi'            };

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_20X_dir{rowTcond})
    
    %Tcell=extractCells();
    Tcell=extractCells(TextraMeasurements);
    
    % add channel-specific information to Tcell
    
    A488_gene=Tcond.A488_gene(rowTcond);
    A555_gene=Tcond.A555_gene(rowTcond);
    A594_gene=Tcond.A594_gene(rowTcond);
    A647_gene=Tcond.A647_gene(rowTcond);
    
    nCell=height(Tcell);
    Tcell=[Tcell,array2table([repmat(A488_gene,nCell,1),repmat(A555_gene,nCell,1),repmat(A594_gene,nCell,1),repmat(A647_gene,nCell,1)],'VariableNames',{'A488_gene','A555_gene','A594_gene','A647_gene'})];
    % save Tcell in extracted directory
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_20X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    
    writetable(Tcell,fullfile(extractedDataDir,'Tcell1.csv')); %for now, save it locally
   
end

%% Extract spot data --> Tspots
% NOTE: need to add the gaussian fitting before this, then also extract those
cd(parentDir)

channelsToExtract={'gfp','gfpa','gfpb','gfpc','tmr','tmra','tmrb','tmrc','alexa','alexaa','alexab','alexac','cy','cya','cyb','cyc'};
intensityCutoffs=repmat(30,1,length(channelsToExtract));

for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_20X_dir{rowTcond})
    
    % get spots
    Tspots=extractSpots(channelsToExtract,'intensitiesAboveGivenCutoff',intensityCutoffs,'getFittedData',false);
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_20X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    
end


%% put all Tcell and Tspot table together into one
condList=1:10;

extractedDataDirForAllCond='paper/extractedData/E157_validation/ReadoutStrip/20X';
if ~isfolder(extractedDataDirForAllCond)
        mkdir(extractedDataDirForAllCond)
end
cd(parentDir)

TcellAll=table();
TspotsAll=table();

for i=1:length(condList)
    condID=condList(i);
    extractedDataDir=Tcond.Extracted_20X_dir{condID};
    cd(extractedDataDir)
    
    Tcell=readtable('Tcell1.csv');
    Tspots=readtable('Tspots1.csv');
    
    TcellAll=[TcellAll;[array2table(repmat(condID,height(Tcell),1),'VariableNames',{'condID'}),Tcell]];
    TspotsAll=[TspotsAll;[array2table(repmat(condID,height(Tspots),1),'VariableNames',{'condID'}),Tspots]];
    
    cd(parentDir)
end

writetable(TcellAll,fullfile(extractedDataDirForAllCond,'Tcell1All.csv'))
writetable(TspotsAll,fullfile(extractedDataDirForAllCond,'Tspots1All.csv'))


%% import data
extractedDataDirForAllCond='paper/extractedData/E157_validation/ReadoutStrip/20X';

TspotsAll=readtable(fullfile(extractedDataDirForAllCond,'Tspots1All.csv'));
TcellAll1=readtable(fullfile(extractedDataDirForAllCond,'Tcell1All.csv'));


% %% check thresholds manually
% cd(parentDir)
% condID=5; % Use pre-strip condIDs
% tiffDir=Tcond.Tiff_20X_dir{condID};
% cd(tiffDir)
% improc2.launchThresholdGUI
% 
% %% Get Tthresh - exposure times to use and threshold to use 
% cd(parentDir)
% %Tplots=readtable(fullfile(extractedDataDir,'Tplots.xlsx')) % has threshold data
% Tthresh=readtable(fullfile(extractedDataDirForAllCond,'Tthresholds.xlsx'));
% temp1=cell(height(Tthresh),1);
% for i=1:height(Tthresh)
%    temp1{i}=str2num(Tthresh.condIDs{i});    
% end
% Tthresh.condIDs=temp1;
% %%% Plot Readout-->Strip data
% Tthresh



%% make TcellAll2,which has numSpots data for before & after stripping
% first 
TcellPlotPaired=table();
condID_readList=Tcond.condID(all([ismember(Tcond.condID,1:10),startsWith(Tcond.readoutCycle,'cycle')],2))';
%
TcellAll2=TcellAll1(:,[2:4,1,5:end]);
nCell=height(TcellAll2)
TcellAll2=[TcellAll2(:,1:end-3),...
                     array2table(cell(nCell,1),'VariableNames',{'A488_channel'}),array2table(nan(nCell,1),'VariableNames',{'A488_thresh'}),array2table(nan(nCell,1),'VariableNames',{'A488_numSpots'}),...
           TcellAll2(:,end-2),...
                     array2table(cell(nCell,1),'VariableNames',{'A555_channel'}),array2table(nan(nCell,1),'VariableNames',{'A555_thresh'}),array2table(nan(nCell,1),'VariableNames',{'A555_numSpots'}),...
           TcellAll2(:,end-1),...
                     array2table(cell(nCell,1),'VariableNames',{'A594_channel'}),array2table(nan(nCell,1),'VariableNames',{'A594_thresh'}),array2table(nan(nCell,1),'VariableNames',{'A594_numSpots'}),...
           TcellAll2(:,end),...
                     array2table(cell(nCell,1),'VariableNames',{'A647_channel'}),array2table(nan(nCell,1),'VariableNames',{'A647_thresh'}),array2table(nan(nCell,1),'VariableNames',{'A647_numSpots'})];
% add channelname and threshold to TcellAll2
for i=1:height(TcellAll2)
    condID=TcellAll2.condID(i);
    arrayNum=TcellAll2.arrayNum(i);
    objNum=TcellAll2.objNum(i);
    
    prefixList={'A488','A555','A594','A647'};
    for ii=1:length(prefixList)
        prefix=prefixList{ii};
        gene=TcellAll2.([prefix,'_gene']){i};
        % A488 gene
        indTthresh=strcmp(Tthresh.gene,gene);
        assert(sum(indTthresh)==1)
        
        if ismember(condID,Tthresh.condIDs{indTthresh})
            channelRajlabimagetools=Tthresh.channelRajlabimagetools(indTthresh);
            threshold=Tthresh.threshold(indTthresh);
            TcellAll2.([prefix,'_channel'])(i)=channelRajlabimagetools;
            TcellAll2.([prefix,'_thresh'])(i)=threshold;
            
            % get numSpots
            TspotsThisCellThisChannel=TspotsAll(all([TspotsAll.condID==condID,TspotsAll.arrayNum==arrayNum,TspotsAll.objNum==objNum,strcmp(TspotsAll.channel,channelRajlabimagetools),TspotsAll.intensities>=threshold],2),:);
            numSpots=height(TspotsThisCellThisChannel);
            TcellAll2.([prefix,'_numSpots'])(i)=numSpots;
        end
    end
    
    
end
%
extractedDataDirForAllCond='extractedData/E157_validation/ReadoutStrip/20X';
writetable(TcellAll2,fullfile(extractedDataDirForAllCond,'TcellAll2.csv'))
%% make TcellPlotPaired, a cell-level data table with paired numSpots data for before & after stripping


for condID_read=condID_readList
    condID_read
    condID_strip=Tcond.condID(Tcond.ReferenceCondID==condID_read)
    
    assert(length(condID_strip)==1)
        
    
    % Find the cells associated with this plotGroupID
    read_indTcell=find(TcellAll2.condID==condID_read);
    strip_indTcell=find(TcellAll2.condID==condID_strip);
    
    %% QC checks
    % number of cells is same?
    if length(read_indTcell)~=length(strip_indTcell)
       error('must have equal number of before and after strip cells') 
    end
    
    % do arrayNum and objNum match
    QC_vect=all([[TcellAll2.arrayNum(read_indTcell),TcellAll2.objNum(read_indTcell),TcellAll2.imageFileSuffix(read_indTcell)]==...
            [TcellAll2.arrayNum(strip_indTcell),TcellAll2.objNum(strip_indTcell),TcellAll2.imageFileSuffix(strip_indTcell)]],2);
    if any(~QC_vect)
        error('the arrayNum, objNum, and imageFileSuffix of the before and after strip cells do not match up')
    end
    % does cellArea match
        % note: if segmentation is too close to edge of image and the XY shift
        % causes it to hit the edge, then the copied-over segmentations will be
        % made small 
    
    if ~isequal(TcellAll2.areaCell(read_indTcell),TcellAll2.areaCell(strip_indTcell))
        error('cell areas do not match')
    end
    
    % gene information QC
    if ~isequal(TcellAll2(read_indTcell,{'A488_gene' 'A555_gene' 'A594_gene' 'A647_gene'}),TcellAll2(strip_indTcell,{'A488_gene' 'A555_gene' 'A594_gene' 'A647_gene'}))
        error('gene labels dont match')
    end
    
    
    %% merge read and strip (r_ and s_) data
    
    vars_Constant=ismember(TcellAll2.Properties.VariableNames,{'arrayNum','objNum','imageFileSuffix','A488_gene' 'A555_gene' 'A594_gene' 'A647_gene'});
    Tcell_read_temp=TcellAll2(read_indTcell,:);   
    Tcell_strip_temp=TcellAll2(strip_indTcell,:);  
    
    % rename most variables to have r_ (readout) or s_ (after stripping)
    Tcell_read_temp.Properties.VariableNames(~vars_Constant)=strcat('r_',TcellAll2.Properties.VariableNames(~vars_Constant));
    Tcell_strip_temp.Properties.VariableNames(~vars_Constant)=strcat('s_',TcellAll2.Properties.VariableNames(~vars_Constant));
    
    vars_JoinExclude=ismember(Tcell_read_temp.Properties.VariableNames,{'A488_gene' 'A555_gene' 'A594_gene' 'A647_gene'});
    TcellPlotPaired_thisPlotGroup=join(Tcell_read_temp(:,~vars_JoinExclude),Tcell_strip_temp(:,~vars_JoinExclude),'Keys',{'arrayNum','objNum','imageFileSuffix'});
    TcellPlotPaired_thisPlotGroup=[TcellPlotPaired_thisPlotGroup(:,1:3),Tcell_read_temp(:,vars_JoinExclude),TcellPlotPaired_thisPlotGroup(:,4:end)];
    
    TcellPlotPaired=[TcellPlotPaired;TcellPlotPaired_thisPlotGroup];

    %% get numSpots from Tspots and do QC checks on them
    TspotsAll_read=TspotsAll(all([TspotsAll.condID==condID_read,strcmp(TspotsAll.channel,channelRajlabimagetools),TspotsAll.intensities>=threshold],2),:);
    TspotsAll_strip=TspotsAll(all([TspotsAll.condID==condID_strip,strcmp(TspotsAll.channel,channelRajlabimagetools),TspotsAll.intensities>=threshold],2),:);
    
    % QC check: make sure cells' arrayNum and objNum match with the ones in TcellAll2
    assert(all(ismember(unique(TspotsAll_read(TspotsAll_read.condID==condID_read,{'arrayNum','objNum'}),'rows'),TcellPlotPaired(TcellPlotPaired.r_condID==condID_read,{'arrayNum','objNum'}))))
    assert(all(ismember(unique(TspotsAll_strip(TspotsAll_strip.condID==condID_strip,{'arrayNum','objNum'}),'rows'),TcellPlotPaired(TcellPlotPaired.s_condID==condID_strip,{'arrayNum','objNum'}))))
    
    % add numSpots into Tcell_read_temp and Tcell_strip_temp
    
    
    
    
end
% reorder
startCols=7;
nMergeCols=width(TcellPlotPaired)-startCols;
newOrder=[1:startCols,reshape([startCols+1:startCols+nMergeCols/2; startCols+1+nMergeCols/2:startCols+nMergeCols],1,nMergeCols)];
assert(length(unique(newOrder))==width(TcellPlotPaired))
TcellPlotPaired=TcellPlotPaired(:,newOrder);
TcellPlotPaired
%

writetable(TcellPlotPaired,fullfile(extractedDataDirForAllCond,'TcellPlotPaired.csv'))


