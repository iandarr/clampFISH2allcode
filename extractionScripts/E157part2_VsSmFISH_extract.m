% E157part2_VsSmFISH_extract
% clampFISH 2.0 vs. smFISH correlation

%% Table of all conditions

%parentDir='/Volumes/IAND_04';
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E157_validation_Tcond.xlsx')
if ~all(Tcond.condID==[1:height(Tcond)]')
    error('condID must be equivalent its row number in this code')
end

condList=[12 14 15 18 20];
%% Segmentation
%navigate to a subfolder and run this. Segment cells.
% do this for 60X
improc2.segmentGUI.SegmentGUI
%% save all segmentations (of source 60X and destination 10X/20X magnifications) to a subfolder for safekeeping
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
%condList=[14];
magnificationList={'60X'};

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
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
end


%%  Process image objects for 60X segmentation data
%condList=[14];
magnificationList={'60X'}; % 60X is processed earlier

sigma=0.5; % use 0.5 for 60X
numLevels=3;


cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        tools = improc2.launchImageObjectTools();
        
        % check if you have saved the segmentations-only .mat files
        if ~isfolder('dataFiles1_segmentations'), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
        
        %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
        improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
        cd(parentDir)
    end
end

%% exclude the bottom slice(s)
%condList=[14];
numBottomSlicesToExclude=1;
magnificationList={'60X'};

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        tools = improc2.launchImageObjectTools();
        
%         % check if you have saved the segmentations-only .mat files
%         if ~isfolder('dataFiles2_processed'), error("save the .mat files with processing done in a subfolder before you begin next steps"),end
        excludeSlices(numBottomSlicesToExclude)
        cd(parentDir)
    end
end
%% save all processed data files for safekeeping
%condList=[14];
magnificationList={'60X'};

subfolderName='dataFiles2_processed';
overwriteFiles=true;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
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
end

%% Thresholding (navigate to data folder, do this for 60X data)
improc2.launchThresholdGUI

%% save thresholded data files for safekeeping
%condList=[14];
magnificationList={'60X'};

subfolderName='dataFiles3_thresholded';
overwriteFiles=true;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
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
end


%% Make Tcell2MapTo20X: Image mapping from source (Ie. 60X) to destination (Ie. 20X and 10X) segmentations
cd(parentDir)
%condList=[14]

sourceMagnification='60X';
sPixelSizeMicrons=2*6.5/60;
sImgSize=[1022 1024];

%destinationMagnification='20X';
destinationMagnificationList={'10X','20X'};
%destinationMagnificationList={'20X'};
%dPixelSizeMicrons=1*6.5/20;
dImgSize=[2044 2048];

RowIsDim='+y';
ColIsDim='-x';

for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(destinationMagnificationList)
        destinationMagnification=destinationMagnificationList{iMag};
        
        
        sourceDir=Tcond.(['Tiff_',sourceMagnification,'_dir']){condID};
        destinationDir=Tcond.(['Tiff_',destinationMagnification,'_dir']){condID};
        
        sourceExtractedDataDir=Tcond.(['Extracted_',sourceMagnification,'_dir']){condID};
        % get destination pixel size
        switch destinationMagnification
            case '20X'
                dPixelSizeMicrons=1*6.5/20;
            case '10X'
                dPixelSizeMicrons=1*6.5/10;
            otherwise
                error('destinationMagnification should be 10X or 20X')
        end
        
        % get a table of Tcell from source folder
        cd(sourceDir)
        Tcell1=extractCells('returnMasks');
        
        cd(parentDir)
        if ~isfolder(sourceExtractedDataDir)
            mkdir(sourceExtractedDataDir)
        end
        Tcell1NoMasks=Tcell1(:,~ismember(Tcell1.Properties.VariableNames,{'cellMaskCropped'}));
        writetable(Tcell1NoMasks,fullfile(sourceExtractedDataDir,filesep,'Tcell1.csv'));
        
        
        % inputs to mapSegsToNewImgs
        % s= source of segmentations (Ie. 60X images with manual segmentations)
        % d= destination of segmentations (Ie. 20X or 10X images)
        
        sImageFileSuffix=Tcell1.imageFileSuffix;
        sBoundingBoxes=[Tcell1.bbox1,Tcell1.bbox2,Tcell1.bbox3,Tcell1.bbox4];
        sMasksCropped=Tcell1.cellMaskCropped;
        s_tempXYdata=readtable(fullfile([sourceDir,filesep,'XYdata.csv']));
        sXimgs=[s_tempXYdata.Xposition];
        sYimgs=[s_tempXYdata.Yposition];
        d_tempXYdata=readtable(fullfile([destinationDir,filesep,'XYdata.csv']));
        dXimgs=[d_tempXYdata.Xposition];
        dYimgs=[d_tempXYdata.Yposition];
        
        XoffsetMicron=Tcond.(['Offset_',sourceMagnification,'to',destinationMagnification,'_x'])(condID);
        YoffsetMicron=Tcond.(['Offset_',sourceMagnification,'to',destinationMagnification,'_y'])(condID);
        
        XoffsetsMicron=repmat(XoffsetMicron,height(Tcell1),1); %replace
        YoffsetsMicron=repmat(YoffsetMicron,height(Tcell1),1);
        
        
        
        % call function to map source segmentations onto destination files
        [dImageFileSuffix,dBoundingBoxes,dMasksCropped,dMasks]=mapSegsToNewImgs(...
            sImageFileSuffix,sBoundingBoxes,sMasksCropped,sXimgs,sYimgs,dXimgs,dYimgs,sImgSize,dImgSize,sPixelSizeMicrons,dPixelSizeMicrons,XoffsetsMicron,YoffsetsMicron,RowIsDim,ColIsDim);
        
        
        % write data00X.mat files in destination folder based on Tcell2
        indSegOut=~isnan(dImageFileSuffix);
        [dArrayNum,dObjNum]=makeDataFilesWithSegmentations(destinationDir,dImageFileSuffix,dMasks);
        
        
        % Make combined table from all magnifications (Tcell2)
        TcellAdd= [...
            array2table(dImageFileSuffix,'VariableNames',{['imageFileSuffix_',destinationMagnification]}),...
            array2table(dArrayNum,'VariableNames',{['arrayNum_',destinationMagnification]}),...
            array2table(dObjNum,'VariableNames',{['objNum_',destinationMagnification]}),...
            array2table(dBoundingBoxes,'VariableNames',{['bbox1_',destinationMagnification],['bbox2_',destinationMagnification],['bbox3_',destinationMagnification],['bbox4_',destinationMagnification]}),...
            array2table(dMasksCropped,'VariableNames',{['masksCropped_',destinationMagnification]}),...
            array2table(dMasks,'VariableNames',{['masks_',destinationMagnification]})];
        
        if iMag==1 %initialize Tcell2
            Tcell2=Tcell1;
            Tcell2.Properties.VariableNames=replace(join([Tcell2.Properties.VariableNames;repmat({'_60X'},1,length(Tcell2.Properties.VariableNames))],1),' ','');
        end
        Tcell2=[Tcell2,TcellAdd];
    end
    
    extractedDataAllMagDir=Tcond.Extracted_AllMagnification_dir{condID};
    if ~isfolder(extractedDataAllMagDir)
        mkdir(extractedDataAllMagDir)
    end
    
    % output Tcell2 but without the masks
    TcellMapNoMasks=Tcell2(:,~ismember(Tcell2.Properties.VariableNames,{'cellMaskCropped_60X','masksCropped_10X','masks_10X','masksCropped_20X','masks_20X'}));
    writetable(TcellMapNoMasks,fullfile(extractedDataAllMagDir,filesep,'TcellMapBetweenMagnifications.csv'))
end

%% pick a row to manually compare the source to destination alignment
% 
% rowsWithDestinationSegmentation=find(~isnan(dArrayNum));
% rowToCompare=rowsWithDestinationSegmentation(1); % 1 for first or end for last
% sourceFileSuffix=sImageFileSuffix(rowToCompare);
% destinationFileSuffix=dArrayNum(rowToCompare);
% % SOURCE
% cd(parentDir)
% cd(sourceDir)
% improc2.segmentGUI.SegmentGUI(sourceFileSuffix)
% % DESTINATION
% cd(parentDir)
% cd(destinationDir)
% improc2.segmentGUI.SegmentGUI(destinationFileSuffix)
% 
% Tcell2
% TcellAdd



%% save all segmentations of destination 10X/20X magnifications to a subfolder for safekeeping
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
%condList=[14];
magnificationList={'10X','20X'};

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
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
end


%%  Process image objects for segmentation destination data (10X and 20X)
%condList=[14];
magnificationList={'10X','20X'}; % 60X is processed earlier

sigma=0.4; % use 0.4 for 10X (0.45NA, 1x1 binning) and 20X (0.75NA, 1x1 binning)
numLevels=3;


cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        tools = improc2.launchImageObjectTools();
        
        % check if you have saved the segmentations-only .mat files
        if ~isfolder('dataFiles1_segmentations'), error("save the .mat files with only segmentations done in a subfolder before you begin processing"),end
        
        %improc2.processImageObjects(); % default sigma=0.5, numLevels=3
        
        improc2.processImageObjects('filterParams', struct('sigma',sigma,'numLevels',numLevels));
        
        cd(parentDir)
    end
end
%% save processed data of destination magnifications for safekeeping
%condList=[14];
magnificationList={'10X','20X'};

subfolderName='dataFiles2_processed';
overwriteFiles=true;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
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
end

%% Extract cell data for 10X and 20X --> Tcell1
% left out extra measurements 
%condList=[14];
magnificationList={'10X','20X'};

% TextraMeasurements=table('Size',[4,4],'VariableNames',{'columnName','area','valueToCompute','channel'},'VariableTypes',{'string','string','string','string'});
% TextraMeasurements.columnName=          {  'avgCellYFPc';      'avgNoncellYFPc'; 'avgNuclearDAPI';  'avgCytoplasmDAPI'};
% TextraMeasurements.area=                {  'cell';             'nonCell';        'nuclear'       ;  'cytoplasmic'      };
% TextraMeasurements.valueToCompute=      {  'average';          'average';        'average'       ;  'average'         };
% TextraMeasurements.channel=             {  'gfpc';             'gfpc';           'dapi'          ;  'dapi'            };

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        Tcell=extractCells();
        %Tcell=extractCells(TextraMeasurements);
        
        % add channel-specific information to Tcell
        %A488_gene=Tcond.A488_gene(condID);
        %A555_gene=Tcond.A555_gene(condID);
        %A594_gene=Tcond.A594_gene(condID);
        %A647_gene=Tcond.A647_gene(condID);
        
        nCell=height(Tcell);
        %Tcell=[Tcell,array2table([repmat(A488_gene,nCell,1),repmat(A555_gene,nCell,1),repmat(A594_gene,nCell,1),repmat(A647_gene,nCell,1)],'VariableNames',{'A488_gene','A555_gene','A594_gene','A647_gene'})];
        % save Tcell in extracted directory
        cd(parentDir)
        
        extractedDataDir=Tcond.(['Extracted_',magnification,'_dir']){condID};
        if ~isfolder(extractedDataDir)
            mkdir(extractedDataDir)
        end
        
        writetable(Tcell,fullfile(extractedDataDir,'Tcell1.csv')); %for now, save it locally
    end
end

%% Extract spot data --> Tspots
%condList=[14];
magnificationList={'60X','10X','20X'};

%channelsToExtract={'gfp','gfpa','gfpb','gfpc','tmr','tmra','tmrb','tmrc','alexa','alexaa','alexab','alexac','cy','cya','cyb','cyc'};
intensityCutoff10X=20;
intensityCutoff20X=30;
intensityCutoff60X=300;
%intensityCutoffs=repmat(intensityCutoff,1,length(channelsToExtract));

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        % get spots
        % Tspots=extractSpots(channelsToExtract,'intensitiesAboveGivenCutoff',intensityCutoffs,'getFittedData',false);
        switch magnification
            case '10X'
                Tspots=extractSpots('all','intensitiesAboveGivenCutoff',intensityCutoff10X,'getFittedData',false);
            case '20X'
                Tspots=extractSpots('all','intensitiesAboveGivenCutoff',intensityCutoff20X,'getFittedData',false);                
            case '60X'
                Tspots=extractSpots('all','intensitiesAboveGivenCutoff',intensityCutoff60X,'getFittedData',false);
        end
            cd(parentDir)
        extractedDataDir=Tcond.(['Extracted_',magnification,'_dir']){condID};
        if ~isfolder(extractedDataDir)
            mkdir(extractedDataDir)
        end
        writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    end
end



%% put all Tcell and Tspot table together into one, from all condID's
% condList=[12 14 15 18 20];
% condList=[14]
magnificationList={'60X','10X','20X'};



extractedDataDirForAllCond='paper/extractedData/E157_validation/VsSmFISH/';
if ~isfolder(extractedDataDirForAllCond)
        mkdir(extractedDataDirForAllCond)
end

cd(parentDir)

TcellAll=table();
TspotsAll=table();
TcellMapAll=table();

for i=1:length(condList)
    condID=condList(i);
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        magnificationNumeric=str2num(replace(magnification,'X',''));
        
        extractedDataDir=Tcond.(['Extracted_',magnification,'_dir']){condID};
        cd(extractedDataDir)
        
        Tcell=readtable('Tcell1.csv');
        Tspots=readtable('Tspots1.csv');
        
        TcellAll=[TcellAll;[array2table([repmat(condID,height(Tcell),1),repmat(magnificationNumeric,height(Tcell),1)],'VariableNames',{'condID','magnification'}),Tcell]];
        TspotsAll=[TspotsAll;[array2table([repmat(condID,height(Tspots),1),repmat(magnificationNumeric,height(Tspots),1)],'VariableNames',{'condID','magnification'}),Tspots]];
        
        cd(parentDir)
    end
    
    % get magnification mapping as well
    TcellMap=readtable(fullfile([Tcond.Extracted_AllMagnification_dir{condID},filesep,'TcellMapBetweenMagnifications.csv']));
    TcellMapAll=[TcellMapAll;[array2table(repmat(condID,height(TcellMap),1),'VariableNames',{'condID'}),TcellMap]];
    
end

writetable(TcellAll,fullfile(extractedDataDirForAllCond,'Tcell1All.csv')) %this doesn't really do anything
writetable(TspotsAll,fullfile(extractedDataDirForAllCond,'Tspots1All.csv'))
writetable(TcellMapAll,fullfile(extractedDataDirForAllCond,'TcellMapBetweenMagnifications.csv'))

%% import data and build TcellMap2 (where thresholds used from Tcond)
extractedDataDirForAllCond='paper/extractedData/E157_validation/VsSmFISH/';

TspotsAll=readtable(fullfile(extractedDataDirForAllCond,'Tspots1All.csv'));
% TcellAll1=readtable(fullfile(extractedDataDirForAllCond,'Tcell1All.csv')); % can get rid of this since it doesn't offer anything
TcellMap1=readtable(fullfile([extractedDataDirForAllCond,filesep,'TcellMapBetweenMagnifications.csv']));

%%% Table of all conditions again
% parentDir='/Volumes/IAND_04';
% %parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
% cd(parentDir)
% % Read in table of conditions (Tcond)
% Tcond=readtable('data/E157_validation/E157_validation.xlsx')
% if ~all(Tcond.condID==[1:height(Tcond)]')
%     error('condID must be equivalent its row number in this code')
% end

%%% build TcellMap2, which has spot counts in it for segmentations that are in common for all magnifications
% this is where thresholds are read in from Tcond

%condList=[14];
%warning('only building table for one condID')

magnificationList={'60X','10X','20X'};
imageFileSuffixVariableNames_AllMags=join([repmat({'imageFileSuffix_'},length(magnificationList),1),magnificationList'],'');
indCellsWithAllMags=~any(isnan(TcellMap1{:,imageFileSuffixVariableNames_AllMags}),2);
TcellMap2=TcellMap1(indCellsWithAllMags,:);

VariableNamesNumSpots=join([repmat({'numSpots_'},length(magnificationList),1),magnificationList'],'')';
VariableNamesChannels=join([repmat({'channel_'},length(magnificationList),1),magnificationList'],'')';
VariableNamesThresholds=join([repmat({'threshold_'},length(magnificationList),1),magnificationList'],'')';
TcellMap2=[TcellMap2,...
    array2table(cell(height(TcellMap2),3),'VariableNames',VariableNamesChannels),...
    array2table(nan(height(TcellMap2),3),'VariableNames',VariableNamesNumSpots),...
    array2table(nan(height(TcellMap2),3),'VariableNames',VariableNamesThresholds)];

cd(parentDir)
for i=1:height(TcellMap2)
    condID=TcellMap2.condID(i);
    
    for iMag=1:length(magnificationList)
        magnification=magnificationList{iMag};
        magnificationNumeric=str2num(replace(magnification,'X',''));
        
        arrayNum=TcellMap2.(['arrayNum_',magnification])(i);
        objNum=TcellMap2.(['objNum_',magnification])(i);
        imageFileSuffix=TcellMap2.(['imageFileSuffix_',magnification])(i);
        
        channel=Tcond.(['spots_',magnification,'_channel']){condID};
        
        % get spots related to this cell, magnidication, and channel using its specific
        % condID, arrayNum, objNum, and channel
        TspotsThisCellMagChannel=TspotsAll(all([TspotsAll.condID==condID,TspotsAll.magnification==magnificationNumeric,TspotsAll.arrayNum==arrayNum,TspotsAll.objNum==objNum,strcmp(TspotsAll.channel,channel)],2),:);
        
        %
        thresholdInput=Tcond.(['spots_',magnification,'_threshold'])(condID);
        if iscell(thresholdInput)
            thresholdInput=thresholdInput{1};
        end
        
        threshold=nan;
        if ischar(thresholdInput)
            if ~isempty(str2num(thresholdInput))
                threshold=str2num(thresholdInput);
            elseif strcmp(thresholdInput,'manual')
                % threshold mode is manual
                thresholdManual=TspotsThisCellMagChannel.threshold(1);
                assert(all(TspotsThisCellMagChannel.threshold==thresholdManual));
                numSpots=sum(TspotsThisCellMagChannel.intensities>=thresholdManual);
                threshold=thresholdManual;
            else
                error('input threshold in Tcond.spots_%s_threshold cannot be this string %s',magnification,thresholdInput)
            end
        elseif isnan(thresholdInput)
            error('input threshold in Tcond.spots_%s_threshold cannot be nan',magnification)
        elseif isnumeric(thresholdInput)
            threshold=thresholdInput;
        end
            numSpots=sum(TspotsThisCellMagChannel.intensities>=threshold);
            % threshold mode is constant
            %error('could not read threshold input in Tcond')
        % store in TcellMap2
        TcellMap2.(['channel_',magnification]){i}=channel;
        TcellMap2.(['numSpots_',magnification])(i)=numSpots;
        TcellMap2.(['threshold_',magnification])(i)=threshold;
        
    end
end

% add on normalizedDistFromCenterOfImage variable

imageHalfDiagonal=sqrt(2)*2048/2;
imageCenterPixelRow=2044/2+0.5;
imageCenterPixelCol=2048/2+0.5;
magnificationList={'10X','20X'};
for iMag=1:length(magnificationList)
    magnification=magnificationList{iMag};
    centroidRow=TcellMap2.(['bbox2_',magnification]) + TcellMap2.(['bbox4_',magnification])/2;
    centroidCol=TcellMap2.(['bbox1_',magnification]) + TcellMap2.(['bbox3_',magnification])/2;
    pixelsFromCenterOfImage=((centroidRow-imageCenterPixelRow).^2 + (centroidCol-imageCenterPixelCol).^2).^0.5;
    normalizedDistFromCenterOfImage=pixelsFromCenterOfImage/imageHalfDiagonal;
    TcellMap2.(['normalizedDistFromCenterOfImage_',magnification])=normalizedDistFromCenterOfImage;
    
    TcellMap2.(['LowerRightness_',magnification])= ((centroidRow - imageCenterPixelRow) + (centroidCol - imageCenterPixelCol))/(imageCenterPixelRow+imageCenterPixelCol);
    TcellMap2.(['Rightness_',magnification])= (centroidCol - imageCenterPixelCol)/(imageCenterPixelCol);
end
TcellMap2

%
writetable(TcellMap2,fullfile([extractedDataDirForAllCond,filesep,'TcellMapBetweenMagnifications2.csv']));






