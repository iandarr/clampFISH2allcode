% E159part2_scan_raw
% high-throughput scan

%% Load Tcond
%clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
Tcond=readtable(fullfile(parentDir,'experiments/E159part2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'experiments/E159part2_scan_Tscan.xlsx'));
Tscan.condIDs=cellfun(@(x) str2num(x),Tscan.condIDs,'UniformOutput',false);
cd(parentDir)
%scanIDsList=Tscan.scanID;
scanIDsList=[1:7]'

diskVolumeName='/Volumes/IAND_06/';
cellposeCommandsFileAllScansDirectory='/Volumes/IAND_06/E159part2_scan/scan';

%subscanGenerationMode='CornerSubscansOnly';

%% make ScanObject.mat file with all rounds, including the XY positions of each tile
% % still need to fix stitching issue - did I shift images in wrong direction.
for iScanID=1:length(scanIDsList)
    
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    condIDsAll=Tscan.condIDs{rowInTscan};
    numRounds=length(condIDsAll);
    Nd2FilepathList=cell(numRounds,1);
    %channelLabels=repmat({{''}},numRounds,1);
    Tchan=table();
    for thisRound=1:length(Nd2FilepathList)
        
        % make Nd2FilepathList
        rowInTcond=find(Tcond.condID==condIDsAll(thisRound));
        [Nd2_filepath,Nd2_filename]=dirAndFilenametext2filepath(Tcond.Nd2_dir{rowInTcond},Tcond.Nd2_filenametext{rowInTcond});
        Nd2FilepathList{thisRound}=Nd2_filepath;
         % make channelLabels input
        channelLabels=Tcond{rowInTcond, {'A488_gene',   'A555_gene',   'A594_gene', 'A647_gene'}};
        channels=                       {'YFP',         'CY3',         'A594',      'CY5'};
        
        %=~ismember(channelLabels,{''}); channels=channels(idxNotEmpty); channelLabels=channelLabels(idxNotEmpty);
        channels=       [{'DAPI','Brightfield'},channels]; channelLabels=  [{'',''},channelLabels]; % append DAPI and Brightfield
        if ismember(scanID,[7 8])
            channelPrefixes=join([repmat(Tcond.cyclePrefixAlt(rowInTcond),1,length(channels));channels],'_',1);
        else
            channelPrefixes=join([repmat(Tcond.cyclePrefix(rowInTcond)   ,1,length(channels));channels],'_',1);
        end
        Tchan=[Tchan; [array2table(repmat(thisRound,length(channels),1),'VariableNames',{'round'}), array2table([channels',channelLabels',channelPrefixes'],'VariableNames',{'channel','channelLabels','channelPrefixes'})]]; % add to channel labels table
    end
    %Tchan
    
    % make R4 (with YFP for cyto segmentations) the reference round
    refRound=4;
    
    %% call makeScanObject
    
    Scan=makeScanObject(Nd2FilepathList,...
       'channelLabels',Tchan,...
       'referenceRound',refRound);%,...
       %'rowDimIsFlippedVsYDim',true,...
       %'cameraAngleOfReferenceRound',-179.7526);%from wellA1, refRound=1, 8pairs: -179.7526 although was using -179.792 for a while % take out 'cameraAngle' pre-specification later on

    
    Scan.scanName=Tscan.scanName{rowInTscan};
    
    scanDir=Tscan.scanDir{rowInTscan};
    if ~isfolder(scanDir)
        mkdir(scanDir)
    end
	save(fullfile(scanDir,filesep,'ScanObject.mat'),'Scan', '-v7.3')
        
end % end scans loop

%% makeStitches

%StitchSubregionSubsetList=[1 1;3 10];
%warning('subregion subset only');
%[a,b]=ind2sub([10,10],1:100);
%SubregionsAll=[b',a'];

%StitchSubregionSubsetList=SubregionsAll(1:16,:);
%StitchSubregionSubsetList=SubregionsAll(17:end,:);

%clear Scan
for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    fprintf('%s\n',repmat('-',1,200))
    fprintf('%s starting makeStitches for scanID=%i at %s  %s\n',repmat('-',1,80),scanID,char(datetime),repmat('-',1,80))
    
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)

	%SubregionArray=[6 6]; % well B3, with resistant cells, had smaller scan (24x24 instead of 39x39)
    %SubregionArray=[10 10];
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    
    scanDir=Tscan.scanDir{rowInTscan};
    load(fullfile(scanDir,filesep,'ScanObject.mat'))
    
    if ismember(scanID,[7,8])
        tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
    else
        tifNameFormat={'R','round','_','channelNames','_','channelLabels'};
    end
    
    makeStitches(Scan,...
        'SubregionArray',SubregionArray,...
        'tifNameFormat',tifNameFormat,...
        'tiffOutParentDir',scanDir,...
        'ScanBounds','innerBoxIntersect')%,...
    %'StitchSubregionSubsetOnly', StitchSubregionSubsetList);
    
end

%% In advance of cellpose segmentations, call mergeImgstoRGB for all images
cd(parentDir)


% scanIDsList=Tscan.scanID
% scanIDsList=[4,5,6]'
%scanIDsList=[7,8]'
% warning('scanIDsList is short list');


convertTo='uint8';
outColors={'g';'b'};

zplanes=1;

for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    % image names to merge
    if ismember(scanID,[7 8])
        imgsForMergePrefixes={'Rcyto_YFP_2000ms';'Rcyto_DAPI_50ms'}; % cytoplasm, nuclear
    else
        imgsForMergePrefixes={'R4_YFP';'R4_DAPI'}; % cytoplasm, nuclear
    end
    imgInputList=join([imgsForMergePrefixes,repmat({'.tif'},size(imgsForMergePrefixes,1),1)],'',2);
    
    
    if ismember(scanID,[7 8])
        outFilenameMerged='Merged_Rcyto_YFP_DAPI.tif';
    else
        outFilenameMerged='Merged_R4_YFP_R4_DAPI.tif';
    end
    
    % scalingMinMax
    cellType=Tscan.cellType{rowInTscan};
    switch cellType
        case 'naive'
            scalingMinMax={[nan nan];[nan nan]};
        case 'resistant'
            scalingMinMax={{0 80};[nan nan]}; % to 80th percentile. Helps because resistant cells are so big that center area is higher in autofluorescence and cellpose makes cells too small. So this merge is more blown out.
        otherwise
            error('can only handle cellType naive and resistant')
    end
    
    FileObj      = java.io.File(diskVolumeName);
    free_Gb   = FileObj.getFreeSpace / 1e9;
    if free_Gb<40
        error('not much disk space left, going to stop this')
    end
    
    fprintf('%s call mergeImgstoRGB for scanID=%i at %s  %s\n',repmat('-',1,20),scanID,char(datetime),repmat('-',1,20))

    
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    numSubregion=length(subregionDirList);
    for iSubregion=1:numSubregion
        
        
        fprintf('    scanID=%i mergeImgstoRGB for subregion %i of %i at %s\n',scanID,iSubregion,numSubregion,char(datetime))
        
        inputImageFullDir=[subregionParentDir,filesep,subregionDirList{iSubregion}];
        
        %outDir=fullfile(inputImageFullDir,filesep,'mergedImgs');
        outDir=inputImageFullDir;
        
        %tiffFileNamesAll=listFilesWithExtension('.tif','directoryPath',inputImageFullDir);
        
            cd(inputImageFullDir)
            imgMerge=mergeImgsToRGB(imgInputList,outColors,'scalingMinMax',scalingMinMax,'zplanes',zplanes,'convertTo',convertTo);
            cd(parentDir)
            
            % write to file
            outFilePath=fullfile(outDir,outFilenameMerged);
            
            imwrite(imgMerge,outFilePath)
            
    end
    
end

%% Optionally add an all-rounds YFP max merge, to help with identifying autoflourescent junk

% scanIDsList=[1,2,3]'
% warning('scanIDsList is short list');

for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    FileObj      = java.io.File(diskVolumeName);
    free_Gb   = FileObj.getFreeSpace / 1e9;
    if free_Gb<40
        error('not much disk space left, going to stop this')
    end
    
    fprintf('%s add YFP max merge for scanID=%i at %s  %s\n',repmat('-',1,20),scanID,char(datetime),repmat('-',1,20))
    
    SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];
    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
    [~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);
    
    numSubregion=length(subregionDirList);
    subregionAbsoluteDirList=join([repmat({subregionParentDir},numSubregion,1),subregionDirList],'',2);

    if ismember(scanID,[7 8])
        R1_TifName='R1_YFP_100ms_UBC.tif';
        R2_TifName='R2_YFP_100ms_UBC.tif';
        R3_TifName='R3_YFP_100ms_UBC.tif';
    else
        R1_TifName='R1_YFP_UBC.tif';
        R2_TifName='R2_YFP_UBC.tif';
        R3_TifName='R3_YFP_UBC.tif';
    end
    
    for iSubregion=1:numSubregion
        cd(subregionAbsoluteDirList{iSubregion});
        
        fprintf('%    scanID=%i add YFP max merge for subregion %i of %i at %s\n',scanID,iSubregion,numSubregion,char(datetime))
        
        imgAllR=uint16(false(imfinfo(R1_TifName).Height,imfinfo(R1_TifName).Width,3));
        imgAllR(:,:,1)=imread(R1_TifName);
        imgAllR(:,:,2)=imread(R2_TifName);
        imgAllR(:,:,3)=imread(R3_TifName);
        imgAllR=max(imgAllR,[],3);
        imwrite(imgAllR,'R123_maxYFP.tif');
    end
end

%% create a .txt file with the commands to call cellpose for all the subregions. For each scan, file outputs to scanDir from Tscan.
% Before you run the text file, you'll need to install cellpose. I used
% conda to install cellpose, and made a conda environment called cellpose.
% Therefore, in terminal I run this first (not including the $):
%   $ conda activate cellpose
% 
% and then run the .txt file, which will be named something like CellposeBatch_cyto_scanID1_for_100_subregions_13-Sep-2021_T20_06_29.txt and contain lines like this:
% 
%   python -m cellpose --dir     /Volumes/IAND_05/E159part2_scan/scan/scan1_WellA1/SubregionArray10x10/Subregion_1_r1_c1 --img_filter Merged_R4_YFP_R4_DAPI --pretrained_model cyto --chan 2 --chan2 3 --fast_mode --no_npy --save_tif --diameter 90
% 
%   python -m cellpose --dir     /Volumes/IAND_05/E159part2_scan/scan/scan1_WellA1/SubregionArray10x10/Subregion_2_r1_c2 --img_filter Merged_R4_YFP_R4_DAPI --pretrained_model cyto --chan 2 --chan2 3 --fast_mode --no_npy --save_tif --diameter 90
% 
% ... and so on for all subregions
% 
% The text file can be run by cd'ing to its directory, Ie:
%   $ cd 
%
% and then just type the full .txt path into Terminal and hitting enter:
%   $ /Volumes/IAND_05/E159part2_scan/scan/scan1_WellA1/CellposeBatch_cyto_scanID1_for_100_subregions_13-Sep-2021_T20_06_29.txt
% 
% or you can cd to the directory and run it like:
%   $ ./MyFilename.txt


cellposeModel='cyto'; % either 'cyto' or 'nuclei'. Make sure cellposeDiam
cytoRGBindex=2; % 2 = green in RGB image
nucleiRGBindex=3; % 3 = blue in RGB image



switch cellposeModel
    case 'cyto'
        chan_str=['--chan ',num2str(cytoRGBindex)];
        chan2_str=['--chan2 ',num2str(nucleiRGBindex)];
        cellposeDiameterIfDrugNaiveCellType=90; % this is the effective 'diameter' of the cell cytoplasm in pixels. 
        cellposeDiameterIfResistantCellType=350;
    case 'nuclei'
        chan_str=['--chan ',num2str(nucleiRGBindex)];
        chan2_str='';
        cellposeDiameterNuclei=50; % this is the typical 'diameter' of the cell nuclei in pixels.
    otherwise
        error('cellposeModel should be cyto or nuclei')
end

tempDatetime=datetime; tempDatetime.Format='ddMMMyyyy_''T''HH_mm_ss'; 
cellposeCommandsFileAllScans=[cellposeCommandsFileAllScansDirectory,filesep,sprintf('Cellpose_allScans_%s.txt',char(tempDatetime))];
fidAllscans=fopen(cellposeCommandsFileAllScans,'w');


for iScanID=1:length(scanIDsList)
    scanID=scanIDsList(iScanID);
    rowInTscan=find(Tscan.scanID==scanID);
    assert(numel(rowInTscan)==1)
    
    if ismember(scanID,[7 8])
        imgNameOfMerged='Merged_Rcyto_YFP_DAPI';
    else
        imgNameOfMerged='Merged_R4_YFP_R4_DAPI';
    end
    
    scanDir=Tscan.scanDir{rowInTscan};
    
    % get cellpose diameter
    cellType=Tscan.cellType{rowInTscan};
    if strcmp(cellposeModel,'cyto') && strcmp(cellType,'naive')
        cellposeDiameter=cellposeDiameterIfDrugNaiveCellType;
    elseif strcmp(cellposeModel,'cyto') && strcmp(cellType,'resistant')
        cellposeDiameter=cellposeDiameterIfResistantCellType;
    elseif strcmp(cellposeModel,'nuclei')
        cellposeDiameter=cellposeDiameterNuclei;
    else
        error('unknown condition')
    end
    
	SubregionArray=[Tscan.SubregionRows(rowInTscan),Tscan.SubregionColumns(rowInTscan)];
    subregionDirStr=['SubregionArray',num2str(SubregionArray(1)),'x',num2str(SubregionArray(2))];

    subregionParentDir=[scanDir,filesep,subregionDirStr,filesep];
    
    subregionDirList=dir(subregionParentDir); subregionDirList=subregionDirList([subregionDirList.isdir]); subregionDirList=subregionDirList(~startsWith({subregionDirList.name},{'.','..'})); subregionDirList={subregionDirList.name}';
	[~,idxSorted]=sort(str2double(regexp(subregionDirList,'\d+','match','once' )));
    subregionDirList=subregionDirList(idxSorted);

    numSubregion=length(subregionDirList);
    subregionAbsoluteDirList=join([repmat({subregionParentDir},numSubregion,1),subregionDirList],'',2);
    
    dirStrPad=max(cellfun(@(x) length(x),subregionAbsoluteDirList));
    
	tempDatetime=datetime; tempDatetime.Format='ddMMMyyyy_''T''HH_mm_ss'; 
    cellposeCommandsFileName=sprintf('Cellpose_%s_scanID%i_for_%i_subregions_%s.txt',cellposeModel,scanID,numSubregion,char(tempDatetime));
    
    fid=fopen(fullfile(scanDir,cellposeCommandsFileName),'w');

    fprintf(fid,'# first you should activate the conda environment in which cellpose in installed. In my case, the conda environment was named cellpose. So the command from Terminal is:\n');
    fprintf(fid,'#    $ conda activate cellpose\n');
    fprintf(fid,'# and then you can call this txt file by typing into Terminal:\n');
    fprintf(fid,'#    $ %s\n\n\n',fullfile(scanDir,filesep,cellposeCommandsFileName));
    

    for iSubregion=1:numSubregion
        
        subregionAbsoluteDir=subregionAbsoluteDirList{iSubregion};
        subregionCellposeStr=sprintf('python -m cellpose --dir %*s --img_filter %s --pretrained_model %s %s %s --fast_mode --no_npy --save_tif --diameter %i',...
                                              dirStrPad,     subregionAbsoluteDir,   imgNameOfMerged,   cellposeModel, chan_str, chan2_str,               cellposeDiameter);
        % file for each well
        fprintf(fid,subregionCellposeStr);
        fprintf(fid,'\n\n');
        
        % file for all wells
        fprintf(fidAllscans,subregionCellposeStr);
        fprintf(fidAllscans,'\n\n');
    end
end

%% call cellpose from Terminal (on Mac) using the output .txt file
% In Terminal, use these commands:
%
%   (base) batfish:~ iandardani$ conda activate cellpose
%   (cellpose) batfish:~ iandardani$ /Volumes/IAND_06/E159part2_scan/scan/Cellpose_allScans_12Jan2022_T21_03_19.txt
%
%
% If you get a message "permission denied" then you may need to move the .txt file to a directory you have execution priveleges with. If that doesn't work try giving execution priveleges to yourself with the Terminal command:   
% 
%   (cellpose) batfish:~ iandardani$ chmod u+x /Volumes/IAND_06/E159part2_scan/scan/Cellpose_allScans_12Jan2022_T21_03_19.txt
% and repeating the above
%% proceed to E159part2_extract


