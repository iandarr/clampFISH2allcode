% E158rep2_E159part1rep2_CellLineAndTissueCompare_raw.m
% Run E158rep2_Tissue_raw first, this uses the stitched outputs from that.
% export a smaller subregion to be processed for spots
% 
% for CellLineVsTissue spot intensity plot:
%   E158rep2_Tissue scan1 (fresh frozen)
%   E158rep2_Tissue scan4 (FFPE)
%   E159part2rep2_scan (cell lines). Already processed.

ImgRegionSize=[500 500]; % pixels rows, column
outImgSize=[2048,2024]; % dentist2 can't handle 500x500 unfortunately
parentDirForOutputOfCroppedRegion='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/rawData/CellLineTissueComparison';
%parentDirForOutputOfCroppedRegion='/Volumes/IAND_09/rawData/CellLineTissueComparison/';
currentDir=pwd;

%% Tissue, FFPE
exportSubDir='TissueFFPE';
iZ=3;
dir_full_subregion=fullfile(parentDirForOutputOfCroppedRegion,filesep,exportSubDir,filesep);
% first output just 1 subregion with single-plane (middle plane)
    SubregionArray=[1 4];
    scanObjectDir='/Volumes/IAND_09/rawData/E158rep2_tissue/scan/scan3_TumorFFPE_CPLX_WM4505-1';
    load(fullfile(scanObjectDir,filesep,'ScanObject.mat'))
    tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
    ScanBounds='innerBoxIntersect';
    makeStitches(Scan,...
        'SubregionArray',SubregionArray,...
        'tifNameFormat',tifNameFormat,...
        'tiffOutParentDir',dir_full_subregion,...
        'ScanBounds',ScanBounds,...
        'zplanes',iZ,...
        'StitchSubregionSubsetOnly',[1,4])
%  next output smaller region
input_dir_full_subregion=fullfile(dir_full_subregion,filesep,'SubregionArray1x4/Subregion_4_r1_c4');
RowCol=[4500 1200]; % top-left pixel's Row, Column
PixelRegion={[RowCol(1),RowCol(1)+ImgRegionSize(1)-1],[RowCol(2),RowCol(2)+ImgRegionSize(2)-1]};
dirOutput=fullfile(parentDirForOutputOfCroppedRegion,filesep,exportSubDir,filesep);
cd(input_dir_full_subregion)
allTifFiles=listFilesWithExtension('.tif');
for iFile=1:length(allTifFiles)
    % these are already single-plane tifs.
    tifName=allTifFiles{iFile};
    img=imread(tifName,'PixelRegion',PixelRegion);
    %imagesc(img)
    % pad with zero pixels since dentist2 can't handle small images
    imgOut=zeros(outImgSize(1),outImgSize(2),'uint16');
    imgOut(1:size(img,1),1:size(img,2))=img;
    imwrite(imgOut,fullfile(dirOutput,filesep,tifName))
end

%% Tissue, Fresh Frozen - get smaller region
exportSubDir='TissueFreshFrozen';
iZ=3;
dir_full_subregion=fullfile(parentDirForOutputOfCroppedRegion,filesep,exportSubDir,filesep);
%dir_full_subregion='/Volumes/IAND_09/rawData/E158rep2_Tissue/scan/scan1_TumorFreshFrozen_Drug/SubregionArray4x4/Subregion_3_r1_c3/';
% first output just 1 subregion with single-plane (middle plane)
input_dir_full_subregion=fullfile(dir_full_subregion,filesep,'SubregionArray4x4/Subregion_3_r1_c3/');
SubregionArray=[4 4];
scanObjectDir='/Volumes/IAND_09/rawData/E158rep2_Tissue/scan/scan1_TumorFreshFrozen_Drug/';
load(fullfile(scanObjectDir,filesep,'ScanObject.mat'))
tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
ScanBounds='innerBoxIntersect';
makeStitches(Scan,...
    'SubregionArray',SubregionArray,...
    'tifNameFormat',tifNameFormat,...
    'tiffOutParentDir',dir_full_subregion,...
    'ScanBounds',ScanBounds,...
    'zplanes',iZ,...
    'StitchSubregionSubsetOnly',[1,3])
%  next output smaller region
RowCol=[2550 1775]; % top-left pixel's Row, Column
PixelRegion={[RowCol(1),RowCol(1)+ImgRegionSize(1)-1],[RowCol(2),RowCol(2)+ImgRegionSize(2)-1]};
dirOutput=fullfile(parentDirForOutputOfCroppedRegion,filesep,exportSubDir,filesep);
cd(input_dir_full_subregion)
allTifFiles=listFilesWithExtension('.tif');
for iFile=1:length(allTifFiles)
    % these are already single-plane tifs.
    tifName=allTifFiles{iFile};
    img=imread(tifName,'PixelRegion',PixelRegion);
    %imagesc(img)
    % pad with zero pixels since dentist2 can't handle small images
    imgOut=zeros(outImgSize(1),outImgSize(2),'uint16');
    imgOut(1:size(img,1),1:size(img,2))=img;
    imwrite(imgOut,fullfile(dirOutput,filesep,tifName))
end

%% CellLine - get region (1 tile) and export to tif
ND2file='/Volumes/IAND_08/rawData/E159part1rep2_ampclick/s17_20X/R8 20220213_201447_672__Well6_ChannelDAPI,YFP,YFP_1,YFP_2,YFP_3,CY3,CY3_1,CY3_2,CY3_3,A594,A594_1,A594_2,A594_3,CY5,CY5_1,CY5_2,CY5_3,Brightfield_Seq0001.nd2';
exportSubDir='CellLine';
dirOutput=fullfile(parentDirForOutputOfCroppedRegion,filesep,exportSubDir,filesep);
tileID=39; % choose a tile
[mSumm,Tstacks,Tplanes]=getMetadata(ND2file);
nc=mSumm.numChannels;
channelNamesIn=mSumm.channelNames';
idx=~contains(channelNamesIn,'_');
temp1=channelNamesIn;
temp1(idx)=append(temp1(idx),'_');
temp2=split(temp1,'_');
temp3=temp2(:,1);
geneNames=replace(temp3,{'YFP','CY3','A594','CY5','DAPI','Brightfield'},{'_UBC','_ITGA3','_FN1','_MITF','',''});
allTifFiles=join([repmat({'R1_'},nc,1),temp3,repmat({'_'},nc,1),replace(cellstr(num2str(mSumm.channelExposureTimesMs')),' ',''),repmat({'ms'},nc,1),geneNames,repmat({'.tif'},nc,1)],'',2);

% first get stack
cd()
for iFile=1:length(allTifFiles)
    thisChannelNameIn=channelNamesIn{iFile};
    tifName=allTifFiles{iFile};
    img = getPlaneFromNd2file(ND2file, tileID, thisChannelNameIn); % this is single-plane dataset
    imwrite(img,fullfile(dirOutput,filesep,tifName))
end