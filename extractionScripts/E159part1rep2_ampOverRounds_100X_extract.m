% E159part1rep2_ampOverRounds_100X_extract
% to determine the spot sizes
% amplification to rounds 1,2,4,6,8,10
% also, amplification to rounds 2,8 without click reaction 

Tcond=readtable('/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/experiments/E159part1rep2_ampclick_Tcond.xlsx');
parentDir='/Volumes/IAND_08';
parentDirForExtractedData='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
%% Segmentation: Navigate to each subfolder and segment cells individually.
% version of rajlabimagetools must have it's getImageFiles.m function with
% the a channels variable that includes 'dapi','gfp','gfpa','gfpb','tmr',tmra','tmrb','alexa','alexaa','alexab','alexac','cya','cyb','cyc','trans'
% ignore the various warnings like "WARNING: Ignoring tmra001.tif file"
improc2.segmentGUI.SegmentGUI
%% save the segmentations to a subfolder for safekeeping
%condList=[1;2;6]
%condList=[3;4]
condList=[5]
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
%condList=Tcond.condID;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowTcond}) % CHOOSE MAGNIFICATION
 
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
%condList=Tcond.condID
%condList=[3]
condList=[5]
display('processing image objects')
origDir=pwd;
sigma=1.3
% theoretical diffraction-limited spot sigmas (approximation)
% 
% set FWHM of gaussian approximation and FWHM of PSF to be equal.
%   FWHM_gaussian = sigma * 2.355
%   FWHM_PSF = =0.51*lambda/NA
% 
%   ---> sigma = 0.51/2.355 * lambda/NA
% 
%   NA=1.45 (100X Apo objective)
%   
% lambda (wavelength) for the channels. Assume mid-point of emission filter: 
%   1. UBC  clampFISH 2.0 with amplifier set 9  labeled in Atto 488  uses 'YFP' filter set (called 'gfp' by rajlabimagetools), which has the emission filter Chroma HQ535/30m --> assume 535nm light
%   2. MITF clampFISH 2.0 with amplifier set 12 labeled in Atto 647N uses 'CY5' filter set (called 'cy'  by rajlabimagetools), which has the emission filter Chroma HQ667/30m --> assume 667nm light
% 
%  sigma (standard deviation of gaussian approximation) in nanometers:
%  1. 'gfp' channel: sigma = 0.51/2.355 * 535nm /1.45 = 79.9nm
%  2. 'cy'  channel: sigma = 0.51/2.355 * 667nm /1.45 = 99.6nm
% 
% camera pixels: with 6.5um camera pixels and with 100X objective (1x1 binning = no binning) = 65nm pixels in object plane.
%
%  sigma in pixels:
%  1. 'gfp' channel: sigma = 79.9nm /(65nm/pixel) = 1.23 pixels
%  2. 'cy'  channel: sigma = 99.6nm /(65nm/pixel) = 1.53 pixels
% 
% chosen sigma = 1.3 pixels

numLevels=3;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    datetime
    fprintf('processing condID=%i\n',condID);
    
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowTcond})
    
    tools = improc2.launchImageObjectTools();
    
    % check if you have saved the segmentations-only .mat files
    if ~isfolder('dataFiles1_segmentations'), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
    
    %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
    improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
    cd(parentDir)
end

% save the processed data files for safekeeping
subfolderName='dataFiles2_processed';
overwriteFiles=true;
%condList=Tcond.condID
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowTcond})
 
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

%% decide on thresholds manually
% navigate to folders with Tiffs in them
improc2.launchThresholdGUI

%% Set uniform threshold for all channels with a value in the threshold column
condList=[1;2;3;4;5;6];
%condList=[6];
cd(parentDir)
channelNamesAll={'gfp','gfpa', 'gfpb', 'cy', 'cya', 'cyb'};
threshVariableNamesAll=join([channelNamesAll',repmat({'_Thresh'},length(channelNamesAll),1)],'',2)';

for i=1:length(condList)
    condID=condList(i)
    rowInTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowInTcond})
        
    thresholds=Tcond{rowInTcond,threshVariableNamesAll};
    channelNames=replace(threshVariableNamesAll,'_Thresh','');
    idx=~isnan(thresholds);
    thresholds=thresholds(idx);
    channelNames=channelNamesAll(idx);

        if any(idx)
            setThreshold(channelNames,thresholds)
        end
    cd(parentDir)
end

%% save the thresholded data files for safekeeping
subfolderName='dataFiles3_thresh1';
overwriteFiles=false;
%condList=Tcond.condID
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowTcond})
 
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

%% Spot fitting


%channelNames={'gfp','gfpa', 'gfpb', 'cy', 'cya', 'cyb'};
channelNames={'gfpb','cyb'};
condList=[1;2;3;4;5;6];
% halfLengthOfRegionToFit. Since 
% 
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    cd(Tcond.Tiff_100X_dir{Tcond.condID==condID})
    fprintf('------- spot fitting condID=%i at %s',condID,char(datetime))
    
    fitSpots(channelNames)
    
    cd(parentDir)
end

%save spot-fitted data files for safekeeping
%save the processed data files for safekeeping AND bring the unfitted
%files to current directory (overwrite)

fprintf('saving fitted data to subfolder')
myDatetime=datetime
subfolderName=sprintf('dataFiles5_fitted_%s',myDatetime);

overwriteFiles=true;
cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowTcond})
 
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

%% Extract Spots --> Tspots1.csv
condList=[1;2;3;4;5;6];
% this works where thresholds are already set
display('------------------------------extracting spots------------------------------')

cd(parentDir)

channelNamesAll={'gfp','gfpa', 'gfpb', 'cy', 'cya', 'cyb'};
threshVariableNamesAll=join([channelNamesAll',repmat({'_Thresh'},length(channelNamesAll),1)],'',2)';

for i=1:length(condList)
    cd(parentDir) % with data
    datetime
    condID=condList(i)
    rowInTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_100X_dir{rowInTcond})

    % which channels to extract
    thresholds=Tcond{rowInTcond,threshVariableNamesAll};
    idx=~isnan(thresholds);
    channelsToExtract=channelNamesAll(idx);

    if any(idx)


        Tspots=extractSpots(channelsToExtract,'intensitiesAboveThreshold','getFittedData',true);

        cd(parentDirForExtractedData) % with extracted data
        extractedDataDir=Tcond.Extracted_100X_dir{rowInTcond};
        if ~isfolder(extractedDataDir)
            mkdir(extractedDataDir)
        end
        writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    end
    cd(parentDir) % with data
end

% Put each condition's Tspots1 table into a single table --> Tspots1All
condList=[1;2;3;4;5;6];
extractedDataDirAllConditions='paper/extractedData/E159part1rep2_ampclick_100X';
cd(parentDirForExtractedData)
TspotsAll=table();
for i=1:length(condList)
    condID=condList(i);
    rowInTcond=find(Tcond.condID==condID);
    condName=Tcond.condName{rowInTcond};
    AmpRound=Tcond.AmpRound(rowInTcond);
    cd(Tcond.Extracted_100X_dir{rowInTcond})
    
    Tspots=readtable('Tspots1.csv');
    Tspots.condID=repmat(condID,height(Tspots),1);
    Tspots.condName=repmat({condName},height(Tspots),1);
    Tspots.AmpRound=repmat(AmpRound,height(Tspots),1);
    Tspots=movevars(Tspots,{'condID','condName','AmpRound'},'Before',1);

    TspotsAll=[TspotsAll;Tspots];
    
    cd(parentDirForExtractedData)
end
writetable(TspotsAll,fullfile(parentDirForExtractedData,filesep,extractedDataDirAllConditions,'Tspots1All.csv'))