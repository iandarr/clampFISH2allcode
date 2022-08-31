% E158rep2_E159part1rep2_CellLineTissueCompare_extract
% will be used in plot comparing spot intensities between cell lines and Fresh Frozen tissue and FFPE tissue
% modeled off of E159part2rep2_scan_extract
% 
% (Run E158rep2_Tissue_raw first)
% 
% for CellLineVsTissue spot intensity plot:
%   E158rep2_Tissue scan1 (fresh frozen)
%   E158rep2_Tissue scan4 (FFPE)
%   E159part2rep2_scan (cell lines). Already processed.

%%
rawDataParentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/rawData/CellLineTissueComparison/';
%rawDataParentDir='/Volumes/IAND_09/rawData/CellLineTissueComparison/';
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';

extractedDataDir=fullfile(parentDir,'extractedData/CellLineTissueComparison',filesep);

% with sigma=0.4
sigma=0.4;
aTrousMinThreshFactor=1.5; % only output spots into spots.csv that are 1.5-fold lower than the threshold (if one is provided), otherwise 1.5-fold lower than the autothreshold for a given block). Eg. if threshold provided is 45, then every spots 30 or greater will be in spots.csv, although all spots <45 will have valid=false
launchGUI=true; % false to run a continuous loop to process spots for all subregions in a batch. Then, afterwards, turn to true to QC check these. The second time you run launchD2ThresholdGUI there will be a spots.csv table in the folder, and it will take these instead of finding spots again
warning('launchGUI=true')
%% for FFPE, call dentist2 
subDir='TissueFFPE';
rawDataDir=fullfile(rawDataParentDir,filesep,subDir,filesep);
cd(rawDataDir)
preStitchedScanFilelist={...
'R1_CY3_1000ms_NGFR.tif',...
'R2_DAPI_50ms.tif',...
'R2_YFP_500ms_UBC.tif',...
'R2_YFP_2000ms_UBC.tif',...
'R2_CY3_1000ms_ITGA3.tif',...
'R2_A594_2000ms_FN1.tif',...
'R3_CY3_2000ms_WNT5A.tif',...
'R3_CY5_500ms_MITF.tif',...
}
channelTypes={'FISH','dapi','FISH','FISH','FISH','FISH','FISH','FISH'};
thresholds=   [        100      100     100     140     100    100    100  ];
% only manually chose threshold for ITGA3 1000ms (thresh=100);
% and mask regions with autofluoresence
h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);

%% for FreshFrozen, call dentist2 
subDir='TissueFreshFrozen';
rawDataDir=fullfile(rawDataParentDir,filesep,subDir,filesep);
cd(rawDataDir)
preStitchedScanFilelist={...
'R1_CY3_1000ms_NGFR.tif',...
'R2_DAPI_50ms.tif',...
'R2_YFP_500ms_UBC.tif',...
'R2_CY3_1000ms_ITGA3.tif',...
'R2_A594_1000ms_FN1.tif',...
'R2_CY5_500ms_EGFR.tif',...
'R3_CY3_1000ms_WNT5A.tif',...
'R3_CY5_500ms_MITF.tif',...
}
channelTypes={'FISH','dapi','FISH','FISH','FISH','FISH','FISH','FISH'};
thresholds=   [100             100     150     100     100    100    100  ];
% only manually chose threshold for ITGA3 1000ms (thresh=100);
h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);    

%% For CellLine, call dentist2
subDir='CellLine';
rawDataDir=fullfile(rawDataParentDir,filesep,subDir,filesep);
cd(rawDataDir)
preStitchedScanFilelist={...
'R1_DAPI_50ms.tif',...
'R1_YFP_500ms_UBC.tif',...
'R1_CY3_1000ms_ITGA3.tif',...
'R1_A594_1000ms_FN1.tif',...
'R1_CY5_500ms_MITF.tif',...
}
channelTypes={'dapi','FISH','FISH','FISH','FISH'};
thresholds=   [        100     180     100   100];
% only manually chose threshold for ITGA3 1000ms (thresh=150);
h=launchD2ThresholdGUI('preStitchedScanFilelist',preStitchedScanFilelist,'launchGUI',launchGUI,'sigma',sigma,'channelTypes',channelTypes,'thresholds',thresholds,'aTrousMinThreshFactor',aTrousMinThreshFactor);    

%% collect spots with actual threshold, fit them, output tables and export to extractedDataDir
TspotsAll=table();
% FFPE

tifFilename='R2_CY3_1000ms_ITGA3';
subDir='TissueFFPE';

fprintf('getting %s spots and fitting them\n',subDir)
Tspots=readtable(fullfile(rawDataParentDir,filesep,subDir,filesep,'spots.csv'));
Tspots.sampleType=repmat({subDir},height(Tspots),1);
idxKeep=all([Tspots.status==1,ismember(Tspots.channel,tifFilename)],2);
Tspots=Tspots(idxKeep,:);
%Tspots=Tspots(1:5:end,:); warning('subsampling');
img=imread(fullfile(rawDataParentDir,filesep,subDir,filesep,[tifFilename,'.tif'])); % only reads single-plane this way
TspotsFit=fitSpotsFromImgWithXYZ(img,Tspots.x,Tspots.y);
TspotsFit.x=[];TspotsFit.y=[];
TspotsFull=[Tspots,TspotsFit];
TspotsAll=[TspotsAll;TspotsFull];

% FreshFrozen
tifFilename='R2_CY3_1000ms_ITGA3';
subDir='TissueFreshFrozen';

fprintf('getting %s spots and fitting them\n',subDir)
Tspots=readtable(fullfile(rawDataParentDir,filesep,subDir,filesep,'spots.csv'));
Tspots.sampleType=repmat({subDir},height(Tspots),1);
idxKeep=all([Tspots.status==1,ismember(Tspots.channel,tifFilename)],2);
Tspots=Tspots(idxKeep,:);
%Tspots=Tspots(1:5:end,:); warning('subsampling');
img=imread(fullfile(rawDataParentDir,filesep,subDir,filesep,[tifFilename,'.tif'])); % only reads single-plane this way
TspotsFit=fitSpotsFromImgWithXYZ(img,Tspots.x,Tspots.y);
TspotsFit.x=[];TspotsFit.y=[];
TspotsFull=[Tspots,TspotsFit];
TspotsAll=[TspotsAll;TspotsFull];

% CellLine
tifFilename='R1_CY3_1000ms_ITGA3';
subDir='CellLine';

fprintf('getting %s spots and fitting them\n',subDir)
Tspots=readtable(fullfile(rawDataParentDir,filesep,subDir,filesep,'spots.csv'));
Tspots.sampleType=repmat({subDir},height(Tspots),1);
idxKeep=all([Tspots.status==1,ismember(Tspots.channel,tifFilename)],2);
Tspots=Tspots(idxKeep,:);
%Tspots=Tspots(1:5:end,:); warning('subsampling');
img=imread(fullfile(rawDataParentDir,filesep,subDir,filesep,[tifFilename,'.tif'])); % only reads single-plane this way
TspotsFit=fitSpotsFromImgWithXYZ(img,Tspots.x,Tspots.y);
TspotsFit.x=[];TspotsFit.y=[];
TspotsFull=[Tspots,TspotsFit];
TspotsAll=[TspotsAll;TspotsFull];

%
writetable(TspotsAll,fullfile(extractedDataDir,filesep,'Tspots.csv'))



