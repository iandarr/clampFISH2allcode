% E159part2rep2_scan_extract
%
% run this after the raw data processing script E159part2rep2_scan_raw.m
%
% it uses these arjunrajlaboratory repos which must be downloaded and on Matlab's path
%   dentist2
%       https://github.com/arjunrajlaboratory/dentist2
%   rajlabimagetools
%       https://github.com/arjunrajlaboratory/rajlabimagetools

%%
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
Tcond=readtable(fullfile(parentDir,'experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'experiments/E159part2rep2_scan_Tscan.xlsx'));
extractedDataDir=fullfile(parentDir,'extractedData/E159part2rep2_scan',filesep);
cd(parentDir)
scanIDsList=Tscan.scanID;
scanIDsList=[1:2]'; % this file for scanID=1,2. (For scanID=3,4 with stripped data and 4th readout cycle data, run E159part2rep2_scanWithStrip_extract.m)
%subscanGenerationMode='CornerSubscansOnly';

%% call dentist2 with preStitchedScanFilelist
% SubregionArray=[10 10]; % now found from Tscan

launchGUI=false; % false to run a continuous loop to process spots for all subregions in a batch. Then, afterwards, turn to true to QC check these. The second time you run launchD2ThresholdGUI there will be a spots.csv table in the folder, and it will take these instead of finding spots again
%warning('launchGUI=true')

preStitchedScanFilelist={...
    'R1_DAPI.tif',...
    'R1_YFP_UBC.tif',...
    'R1_CY3_WNT5A.tif',...
    'R1_A594_DDX58.tif',...
    'R1_CY5_AXL.tif',...    
    'R2_DAPI.tif',...
    'R2_YFP_UBC.tif',...    
    'R2_CY3_NGFR.tif',...
    'R2_A594_FN1.tif',...
    'R2_CY5_EGFR.tif',...
    'R3_DAPI.tif',...
    'R3_YFP_UBC.tif',...
    'R3_CY3_ITGA3.tif',...
    'R3_A594_MMP1.tif',...
    'R3_CY5_MITF.tif',...
    'R4_DAPI.tif',...
    'R4_YFP.tif',...
    'R4_Brightfield.tif',...
 	'Merged_R4_YFP_R4_DAPI_cp_masks.tif',...
	'R123_maxYFP.tif',...
    };
channelTypes={...
    'other','FISH','FISH','FISH','FISH',...
    'other','FISH','FISH','FISH','FISH',...
    'other','FISH','FISH','FISH','FISH',...
    'dapi','other','other',...
    'other','other'};

% with sigma=0.4
sigma=0.4;
%              R1                      R2                      R3   
%              UBC   WNT5A DDX58 AXL   UBC   NGFR  FN1   EGFR  UBC   ITGA3 MMP1  MITF
thresholds=   [37    90    30    60    40    45    55    45    40    45    45    35]; % modified 21-Sep-2021. Even though AXL thresh doesn't match with the E157 threshold 
aTrousMinThreshFactor=1.5; % only output spots into spots.csv that are 1.5-fold lower than the threshold (if one is provided), otherwise 1.5-fold lower than the autothreshold for a given block). Eg. if threshold provided is 45, then every spots 30 or greater will be in spots.csv, although all spots <45 will have valid=false

cd(parentDir); 
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    fprintf('%s launchD2ThresholdGUI for scanID=%i at %s  %s\n',repmat('-',1,20),scanID,char(datetime),repmat('-',1,20))
        
    numSubregion=length(subregionDirList);
    
for iSubregion=1:numSubregion
    
    subregionDir=subregionDirList{iSubregion};
    cd([subregionParentDir,filesep,subregionDir])
    
    fprintf(' -------- scanID=%i, iSubregion=%i (of %i), subregionDir=%s: running launchD2ThresholdGUI starting %s --------\n',scanID,iSubregion,numSubregion,subregionDir,char(datetime))
    %h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',true,'channelTypes',channelTypes,'thresholds',thresholds)
    h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);
    
    %warning('add back in cd parentDir')
    cd(parentDir); % can only cd here if launchGUI=false otherwise it'll get out of it's current directory and not like savingworks
end
end

%% make spotsValidWithSeg.csv. dentist2 is only able to assign spots to nearest nucleus. Now assign spots to cytoplasm

cd(parentDir);
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    %warning('only looking at some subregions')
    %subregionDirList=subregionDirList(ismember(subregionDirList,'Subregion_45_r5_c5')); % select just this one for now
    
    numSubregion=length(subregionDirList);
    
for iSubregion=1:numSubregion

    subregionDir=subregionDirList{iSubregion};
    fprintf('scanID=%i: assigning spots to cytoplasm for subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))
    cd([subregionParentDir,filesep,subregionDir])
    
        
    opts = detectImportOptions('spots.csv');
    opts = setvartype(opts, {'spotID', 'x', 'y', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
    opts = setvartype(opts, 'channel', 'string');
    opts = setvartype(opts, {'intensity', 'status'}, 'uint16'); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
    
    spotTableRaw=readtable('spots.csv',opts);
    spotTableRaw = convertvars(spotTableRaw ,'status', 'logical');
    spotTableValid=spotTableRaw(spotTableRaw.status,:); % only spots with status=true, which means they are above threshold and not masked
    spotTableValid=convertvars(spotTableValid,'channel','categorical');
    
    % assign spot to segmentation, which are given non-zero uint16 numbers
    segmentationImg=imread('Merged_R4_YFP_R4_DAPI_cp_masks.tif');
    
    % [X, Y] = ind2sub(size(bw),spotInds)
    % IND = sub2ind(SIZ,I,J) % where I=row and J=col
    ind=sub2ind(size(segmentationImg),spotTableValid.x,spotTableValid.y); % weirdly x is row, y is col
    spotTableValid.segID=segmentationImg(ind);
    
    spotTableValidWithSeg=spotTableValid(spotTableValid.segID>0,:);
    writetable(spotTableValidWithSeg,'spotsValidWithSeg.csv');
    cd(parentDir);
end
end

%% output spots above threshold (not used in main analyses)
channelsWithFISH={'R1_YFP_UBC',...
    'R1_CY3_WNT5A',...
    'R1_A594_DDX58',...
    'R1_CY5_AXL',...    
    'R2_YFP_UBC',...    
    'R2_CY3_NGFR',...
    'R2_A594_FN1',...
    'R2_CY5_EGFR',...
    'R3_YFP_UBC',...
    'R3_CY3_ITGA3',...
    'R3_A594_MMP1',...
    'R3_CY5_MITF'};

Tthresh=table();
Tthresh.channel=channelsWithFISH';
Tthresh.channel=string(Tthresh.channel);
Tthresh.threshold=thresholds';


cd(parentDir);
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    %warning('only looking at some subregions')
    %subregionDirList=subregionDirList(ismember(subregionDirList,'Subregion_45_r5_c5')); % select just this one for now
    
    numSubregion=length(subregionDirList);
    
for iSubregion=1:numSubregion

    subregionDir=subregionDirList{iSubregion};
    fprintf('scanID=%i: getting spots to make spotTableAboveThreshIgnoringValidity subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))
    cd([subregionParentDir,filesep,subregionDir])
    
        
    opts = detectImportOptions('spots.csv');
    opts = setvartype(opts, {'spotID', 'x', 'y', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
    opts = setvartype(opts, 'channel', 'string');
    opts = setvartype(opts, {'intensity', 'status'}, 'uint16'); %For some reason, when I set 'status' to 'logical' they all go to false. So doing this instead
    
    spotTableRaw=readtable('spots.csv',opts);
    
    % apply thresholds
    spotTableRaw=join(spotTableRaw,Tthresh,'Keys',{'channel'});
    spotTableRaw.aboveThresh=spotTableRaw.intensity>=spotTableRaw.threshold;
    spotTableAboveThreshIgnoringValidity=spotTableRaw(spotTableRaw.aboveThresh,:);
        
    writetable(spotTableAboveThreshIgnoringValidity,'spotTableAboveThreshIgnoringValidity.csv');
    
    cd(parentDir);
end
end
%% Make Tcell1.csv: Add up the spots in each cell to create Tcell for each subregion

cd(parentDir);
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    
    numSubregion=length(subregionDirList);
    
for iSubregion=1:numSubregion
    
    subregionDir=subregionDirList{iSubregion};
    fprintf('scanID=%i: Adding up spots in each cell to make Tcell for subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))
    cd([subregionParentDir,filesep,subregionDir])
    
    opts = detectImportOptions('spotsValidWithSeg.csv');
    opts = setvartype(opts, {'spotID', 'x', 'y', 'nearestNucID', 'maskID', 'distanceToNuc'}, 'single');
    opts = setvartype(opts, {'segID'}, 'uint32');
    opts = setvartype(opts, 'channel', 'categorical');
    opts = setvartype(opts, {'intensity', 'status'}, 'uint16');
    
    Tspots=readtable('spotsValidWithSeg.csv',opts);
    TspotsSumm = groupsummary(Tspots, {'channel', 'segID'}, 'IncludeEmptyGroups', true);
    %channelNames=categories(TspotsSumm.channel);
    TcellWithSpots=unstack(TspotsSumm,'GroupCount','channel');
    
    % if there are no spots for a given channel, include these channels in the table 
    fishChannels={'R1_YFP_UBC','R1_CY3_WNT5A','R1_A594_DDX58','R1_CY5_AXL','R2_YFP_UBC','R2_CY3_NGFR','R2_A594_FN1','R2_CY5_EGFR','R3_YFP_UBC','R3_CY3_ITGA3','R3_A594_MMP1','R3_CY5_MITF'};
    idxFISHToAdd=~ismember(fishChannels,TcellWithSpots.Properties.VariableNames);
    if any(idxFISHToAdd)
        strShow=join(fishChannels(idxFISHToAdd),' ');
        warning('NO FISH SPOTS DETECTED IN %i channels: %s',sum(idxFISHToAdd),strShow{1});
        % pad table with new variables
        Tpad=array2table(zeros(height(TcellWithSpots),sum(idxFISHToAdd),'uint16'),'VariableNames',fishChannels(idxFISHToAdd));
        TcellWithSpots=[TcellWithSpots,Tpad];
    end
    TcellWithSpots=movevars(TcellWithSpots,fishChannels,'After','segID'); % reorder
    
    % if there were some segID with no spots from any channel, need to add
    % it in
    numSegs=max(imread('Merged_R4_YFP_R4_DAPI_cp_masks.tif'),[],'all');
    Tcell=array2table(uint32((1:numSegs)'),'VariableNames',{'segID'});
    Tcell=fillmissing(outerjoin(Tcell,TcellWithSpots,'MergeKeys',true),'constant',0); % since TcellWithSpots is missing some cells, now make Tcell where missing cells have zeros for all channels
    
    fullPathName=fullfile([subregionParentDir,filesep,subregionDir,filesep],'Tcell1.csv');
    writetable(Tcell,fullPathName)
    
    cd(parentDir)
end
    
end
%% Make Tcell2.csv: Add cell properties to Tcell: 'Area', 'Centroid', 'BoundingBox','Circularity','Eccentricity','Orientation','Perimeter','Solidity'

cd(parentDir);
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    UmPerPixel=Tscan.UmPerPixel(rowInTscan);
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    numSubregion=length(subregionDirList);
    
for iSubregion=1:numSubregion
    subregionDir=subregionDirList{iSubregion};
    fprintf('scanID=%i: adding cell area properties to Tcell for subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))
    
    
    maskPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,'Merged_R4_YFP_R4_DAPI_cp_masks.tif');
    maskImg=imread(maskPath);
    stats=regionprops(maskImg,{'Area', 'Centroid', 'BoundingBox','Circularity','Eccentricity','Orientation','Perimeter','Solidity'});
    Tprops=struct2table(stats);
    Tprops=splitvars(Tprops,{'Centroid','BoundingBox'},'NewVariableNames',{{'centCol','centRow'},{'bboxCol','bboxRow','bboxWidth','bboxHeight'}});
    
    outlinesPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,'Merged_R4_YFP_R4_DAPI_cp_outlines.txt');
    Tprops.Area=UmPerPixel^2 * Tprops.Area;
    Tprops=renamevars(Tprops,'Area','Area_Um2');
    Tprops.Perimeter=UmPerPixel * Tprops.Perimeter;
    Tprops=renamevars(Tprops,'Perimeter','Perimeter_Um');
    
    
    inTablePath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'Tcell1.csv');
    Tcell=readtable(inTablePath);
    
    % add Tprops to Tcell
    if height(Tcell)~=height(Tprops)
        fprintf('WARNING: Tprops and Tcell do not have same height. Check this out\n'); % breakpoint
    end
    Tcell=[Tcell,Tprops]; % Tprops must be in same order as Tcell
    
    outTablePath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'Tcell2.csv');
    writetable(Tcell,outTablePath)
    
    
end
    
end


%% combine all the Tcell tables into one

%scanIDsList=Tscan.scanID
%scanIDsList=Tscan.scanID;
%scanIDsList=[1:6]';

%scanIDsList=1;
%warning('scan ID list is smaller subset')

TcellAll=table();

cd(parentDir);
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    scanWell=Tscan.scanWell{rowInTscan};
    
    cellType=Tscan.cellType{rowInTscan};
    assert(any(ismember(cellType,{'naive','resistant'})))
    isResistant=strcmp(cellType,'resistant');
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    numSubregion=length(subregionDirList);
    
    
for iSubregion=1:numSubregion
    subregionDir=subregionDirList{iSubregion};
    fprintf('scanID=%i: adding subregion Tcell to TcellAll from subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))
    cd([subregionParentDir,filesep,subregionDir])
    
    opts = detectImportOptions('Tcell2.csv');
    opts = setvartype(opts, {'segID'}, 'uint32');
    
    Tcell=readtable('Tcell2.csv',opts);
    numCellNew=height(Tcell);
    % append subregion-specific information
    Tcell.scanID=repmat(scanID,numCellNew,1);
    Tcell.scanWell=categorical(repmat({scanWell},numCellNew,1));
    Tcell.isResistant=repmat(isResistant,numCellNew,1);
    
    % get subregion from subregionDir name
    temp1=split(subregionDir,'_');
    temp2=str2double(replace(split(temp1(2:end),'_'),{'r','c'},''));
    subregionInd=temp2(1);subregionRow=temp2(2);subregionCol=temp2(3);
    Tcell.subregionInd=repmat(subregionInd,numCellNew,1);
    Tcell.subregionRow=repmat(subregionRow,numCellNew,1);
    Tcell.subregionCol=repmat(subregionCol,numCellNew,1);
    Tcell=movevars(Tcell,{'scanID','scanWell','subregionInd','subregionRow','subregionCol'},'Before','segID');
    TcellAll=[TcellAll;Tcell];
    cd(parentDir)
end

end
% add unique cellIDs
TcellAll.cellID=[1:height(TcellAll)]';
TcellAll=movevars(TcellAll,{'cellID'},'Before','scanID');
%TcellAll=addvars(TcellAll,[1:height(TcellAll)]','NewVariableNames','cellID','Before',1);

% add uniqueSubregionInd (scanWell and subregionInd define this)
Ttemp=unique(TcellAll(:,{'scanWell','subregionInd'}),'rows');
Ttemp.uniqueSubregionInd=[1:height(Ttemp)]';
TcellAll=join(TcellAll,Ttemp,'Keys',{'scanWell','subregionInd'});% add uniqueSubregionInd 
TcellAll=movevars(TcellAll,'uniqueSubregionInd','After','subregionCol');

fullPathName=fullfile(extractedDataDir,'TcellAll1.csv');
writetable(TcellAll,fullPathName)


%% Apply QC filters TcellAll1.csv --> TcellAll1_withQC.csv and TcellAll2.csv
TcellAll=readtable(fullfile(extractedDataDir,'TcellAll1.csv'));
%TcellAll.isResistant=logical(TcellAll.isResistant);
%TcellAll=TcellAll(TcellAll.isResistant,:);
%warning('remove ~isResistant filter')

% QC critera 1: only subregions with meanUBC > = subregionMinMeanUBC
Tsr=groupsummary(TcellAll(:,{'uniqueSubregionInd','isResistant','subregionInd','subregionRow','subregionCol','scanID','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC'}),'uniqueSubregionInd','mean');
Tsr=renamevars(Tsr,{'mean_subregionInd','mean_isResistant','mean_subregionRow','mean_subregionCol','mean_scanID','GroupCount'},{'subregionInd','isResistant','subregionRow','subregionCol','scanID','numCellsBeforeQC'});
Tsr=join(Tsr,Tscan(:,{'scanID','scanWell'}),'Keys',{'scanID'});
Tsr=movevars(Tsr,'scanWell','After','uniqueSubregionInd');
subregionMinMeanUBC=25; % min of the average of all UBC per cell for each subregion. If you masked a lot of spots (and their area is now zero) then this would produce a lot of zero-count cells, and would bring down this average
% ^ would probably need to be increased for well B3 (scan 6), but this was imaged in a smaller region in the middle of the well, so coverslip edge  problem not applicable

Tsr.passSubregion=all([Tsr.mean_R1_YFP_UBC,Tsr.mean_R2_YFP_UBC,Tsr.mean_R3_YFP_UBC]>=subregionMinMeanUBC,2);

TcellAll=join(TcellAll,Tsr(:,{'uniqueSubregionInd','passSubregion'}));

% QC criteria 2: For all 3 rounds, UBC count >= minUBC
%   Note: cells can have zero UBC if the spots were masked in dentist2, or
%   if truly no spots above threshold were in that segmentation

UBC3Rounds=TcellAll{:,{'R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC'}};

minUBC=4; % per cell;
passMinUBC=all(UBC3Rounds>=minUBC,2);
TcellAll.passMinUBC=passMinUBC;

% QC criteria 3: UBC >= minUBCPerUm2 * CellArea
%   Note: larger cells have more UBC. If it doesn't, it's likely out of focus
minUBCPerUm2=0.025;

passMinUBCPerArea = min(UBC3Rounds,[],2)>= minUBCPerUm2 * TcellAll.Area_Um2;
TcellAll.passMinUBCPerArea=passMinUBCPerArea;

% QC critera 4: Percent UBC error cutoff
minPcrErrUBC=50; % 50 = UBC from all rounds should have less than or equal to 50% difference from median of the 3 rounds
passPctUBC=all((100 * abs(UBC3Rounds./median(UBC3Rounds,2) - 1)) <= minPcrErrUBC,2);
TcellAll.passPctUBC=passPctUBC;

% passes all QC checks
passAllQC=all([TcellAll.passSubregion, TcellAll.passMinUBC, TcellAll.passMinUBCPerArea, TcellAll.passPctUBC],2);
TcellAll.passAllQC=passAllQC;

% add number of cells passe to Tsr
temp=groupsummary(TcellAll(:,{'uniqueSubregionInd','passAllQC'}),'uniqueSubregionInd','sum');
temp=renamevars(temp,'sum_passAllQC','numCellsAfterQC');
Tsr=join(Tsr,temp(:,{'uniqueSubregionInd','numCellsAfterQC'}));
Tsr.pctPassingQC=100*Tsr.numCellsAfterQC./Tsr.numCellsBeforeQC;
Tsr=movevars(Tsr,{'numCellsAfterQC','pctPassingQC'},'After','numCellsBeforeQC');

% statistics for each cell type
for i=1:2
    if i==1
        cellType='naive';
        TcellAll_subset=TcellAll(TcellAll.isResistant==0,:);
        Tsr_subset=Tsr(Tsr.isResistant==0,:);
    elseif i==2
        cellType='resistant';
        TcellAll_subset=TcellAll(TcellAll.isResistant==1,:);
        Tsr_subset=Tsr(Tsr.isResistant==1,:);
    end
    
numOrig=height(TcellAll_subset);
numPassedAllQC=sum(TcellAll_subset.passAllQC);
numPassedSubregion=sum(TcellAll_subset.passSubregion);
uniqueSubregionInds=unique(TcellAll_subset.uniqueSubregionInd);
numGoodSubregion=sum(Tsr_subset.passSubregion);
numTotSubregion=height(Tsr_subset);

fprintf('------------- statistics for cellType=%s ---------------\n',cellType)
fprintf(['QC step 1 (subregion-level): with subregionMinMeanUBC=%.2f, %i subregions (of %i) passed (%i of %i cells, or %.2f%%)\n',...
         'QC step 2 (cell-level):      %i of %i cells (%.2f%%) passed cell-specific QC checks \n',...
         'Overall: %i (%.2f%%) of cells passed all QC\n'],...
                           subregionMinMeanUBC, numGoodSubregion,numTotSubregion, numPassedSubregion, numOrig, 100*numPassedSubregion/numOrig,...
                           numPassedAllQC,numPassedSubregion,100*numPassedAllQC/numPassedSubregion,...
                           numPassedAllQC,100*numPassedAllQC/numOrig)

end
%groupsummary(TcellAll(:,{'subregionInd','R1_YFP_UBC'}),'subregionInd','mean')
%groupsummary(TcellAll(:,{'subregionInd','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC','passMinUBC','passPctUBC','passAllQC'}),'subregionInd','mean')
%groupsummary(TcellAll(:,{'uniqueSubregionInd','scanID','subregionInd','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC','passMinUBC','passPctUBC','passAllQC'}),'uniqueSubregionInd','mean')

writetable(Tsr,fullfile(extractedDataDir,'Tsubregions.csv'))
writetable(TcellAll,fullfile(extractedDataDir,'TcellAll1_withQC.csv'))

TcellAll2=TcellAll(TcellAll.passAllQC,:); % just what passes QC
writetable(TcellAll2,fullfile(extractedDataDir,'TcellAll2.csv'))
