% E159part1rep2_ampOverRounds_extract
% amplification to rounds 1,2,4,6,8,10
% also, amplification to rounds 2,8 without click reaction


% two parts to this experiment, with corresponding analyses:
% E159part1 = fold-change/click in 8-well chambers
% E159part2 = high-throughput cell profiling in 6-well chambers
% 

Tcond=readtable('/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/experiments/E159part1rep2_ampclick_Tcond.xlsx');
parentDir='/Volumes/IAND_08';

%% In advance of cellpose segmentations, call mergeImgstoRGB for all images
cd(parentDir)
%inputDirList=Tcond.Tiff_60X_dir;
inputDirList=Tcond.Tiff_60X_dir(ismember(Tcond.condID,[1,2,3,4,5,6,7,8]));

%inputDirList=  {'E159_scan_amp/s17_8well_60X_stacks/c1_R1_20210612_173333_896__Well2';...
%                'E159_scan_amp/s17_8well_60X_stacks/c3_R10_20210612_192114_033__Well3'};


imgInputPrefixes={'gfpe';'dapi'};

convertTo='uint8';
outColors={'g';'b'};
scalingMinMax={{0 99.5};[1000 30000]};
zplanes=4:8;

for iInputDir=1:length(inputDirList)
    
    
    inputImageDir=inputDirList{iInputDir};
    datetime
    inputImageDir
    
    outDir=fullfile(inputImageDir,filesep,'mergedImgs');
    
    tiffFileNamesAll=listFilesWithExtension('.tif','directoryPath',inputImageDir);
    tiffFileNamesCytoplasm=tiffFileNamesAll(startsWith(tiffFileNamesAll,imgInputPrefixes{1}));
    
    for iFile=1:length(tiffFileNamesCytoplasm)
        % get names of cytoplasm (yfp) and nuclear (dapi) channels
        filenameCytoplasm=tiffFileNamesCytoplasm{iFile};
        numericPortionOfFilename=filenameCytoplasm(regexp(filenameCytoplasm,'\d'));
        filenameDAPI=[imgInputPrefixes{2},numericPortionOfFilename,'.tif'];
        imgInputList={filenameCytoplasm;filenameDAPI};
        
        cd(inputImageDir)
        imgMerge=mergeImgsToRGB(imgInputList,outColors,'scalingMinMax',scalingMinMax,'zplanes',zplanes,'convertTo',convertTo);
        cd(parentDir)
        
        % write to file
        outFilePath=fullfile(outDir,filesep,['Merged',numericPortionOfFilename,'.tif']);
        
        % make output directory if needed
        [tempOutfileDir,~,~]=fileparts(outFilePath);
        if ~isfolder(tempOutfileDir)
            mkdir(tempOutfileDir);
        end
        tempOutfileDir=[];
        
        imwrite(imgMerge,outFilePath)
    end
    
end
% %% delete mergedImgs directories
% cd(parentDir)
% inputDirList=Tcond.Tiff_60X_dir;
% 
% for iInputDir=1:length(inputDirList)
% 
%     inputImageDir=inputDirList{iInputDir};
%     outDir=fullfile(inputImageDir,filesep,'mergedImgs');
% 
%     if isdir(outDir)
%     %    rmdir(outDir,'s')
%     fprintf('%s is a folder\n',inputImageDir)
%     else
%         fprintf('%s is NOT a folder OH NO\n',inputImageDir)
%         cd(inputImageDir)
%         display('but I can cd there')
%         
%     end
%     cd(parentDir)
% end
%% Run cellpose on merged (YFP=green + DAPI=blue) RGB images
% Do this part with batch file called from terminal, see google doc
% this will segment a single directory:
%
%
% (base) batfish:~ iandardani$ conda activate cellpose
% (cellpose) batfish:~ iandardani$ /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/batchCallCellpose.txt
%
% where batchCallCellpose.txt is:
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R1_20220213_212821_051__Well4/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R2_20220213_212821_051__Well2/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R4_20220213_212821_051__Well3/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R6_20220213_212821_051__Well8/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R8_20220213_212821_051__Well6/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R10_20220213_212821_051__Well7/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R2_NoCu_20220213_212821_051__Well1/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150
% python -m cellpose --dir /Volumes/IAND_08/rawData/E159part1rep2_ampclick/s18_60X/R8_NoCu_20220213_212821_051__Well5/mergedImgs --pretrained_model cyto --chan 2 --chan2 3 --save_png --fast_mode --diameter 150

%% Convert cellpose segmentation files and turn into rajlabimagetools data files, equivalent to the output of segmentation GUI (Eg. mm --> data001.mat)
cd(parentDir)
%condList=Tcond.condID;
condList=[1:8];
%warning('make condList actually full')

for i=1:length(condList)
    condID=condList(i);
    tiffDir=Tcond.Tiff_60X_dir{Tcond.condID==condID};
    cellposeSegmentationDir=fullfile(tiffDir,filesep,'mergedImgs');
    
    % get list of segmentation png files
    allPngFiles=listFilesWithExtension('.png','directoryPath',cellposeSegmentationDir);
    segmentationPngFiles=allPngFiles(endsWith(allPngFiles,'_cp_masks.png'));
    nSegmentationFiles=length(segmentationPngFiles);
    segmentationFileNums=str2num(str2mat(replace(segmentationPngFiles,{'Merged','_cp_masks.png'},'')));
    %segmentationPngFiles=replace(segmentationPngFiles,'_cp_masks.png','');
    
    % get list of tif files to check that DAPI file numbers match
    % segmentation numbers
    allTiffFiles=listFilesWithExtension('.tif','directoryPath',tiffDir);
    dapiTiffFiles=allTiffFiles(startsWith(lower(allTiffFiles),'dapi'));
    dapiFileNums=str2num(str2mat(replace(dapiTiffFiles,{'dapi','DAPI','.tif'},'')));
    %assert(length(nSegmentationFiles)==length(dapiTiffFiles));
    assert(isequal(segmentationFileNums,dapiFileNums))
    
    
    
    % get cellpose mask png files with each FOV's segmentations. Put masks in cell array
    maskImgBinaryList=cell(0,1);
    arrayNumList=nan(0,1);
    for iImgFiles=1:length(segmentationPngFiles)
        maskFile=fullfile(cellposeSegmentationDir,segmentationPngFiles{iImgFiles});
        maskImgMultipleSegmentations=imread(maskFile);
        numSegmentationsInThisImg=double(max(max(maskImgMultipleSegmentations)));
        
        thisImgMaskImgBinaryList=cell(numSegmentationsInThisImg,1);
        for iSegmentation=1:numSegmentationsInThisImg
            thisImgMaskImgBinaryList{iSegmentation}=maskImgMultipleSegmentations==iSegmentation;
        end
        thisImgArrayRepmat=repmat(segmentationFileNums(iImgFiles),numSegmentationsInThisImg,1);
        
        % build complete list of masks and associated arrayNum
        maskImgBinaryList=[maskImgBinaryList;thisImgMaskImgBinaryList];
        arrayNumList=[arrayNumList;thisImgArrayRepmat];
    end

    [arrayNum,objNum]=makeDataFilesWithSegmentations(tiffDir,arrayNumList,maskImgBinaryList,'ignoreSmallNonconnectedSegmentations',true);
    
    
    
end


%% save the segmentations to a subfolder for safekeeping
condList
subfolderName='dataFiles1_segmentations';
overwriteFiles=true;
%condList=Tcond.condID;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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
condList
display('processing image objects')
origDir=pwd;
sigma=0.5;
numLevels=3;

cd(parentDir)
for i=1:length(condList)
    condID=condList(i);
    datetime
    fprintf('processing condID=%i\n',condID);
    
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
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
    cd(Tcond.Tiff_60X_dir{rowTcond})
 
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
% left out extra measurements since not all have nuclear mask
% extimating total will take 32min. 12s/img * 8cond * 20imgs/cond
%condList=Tcond.condID
%condList=[1, 3, 5]
%condList=[2]
%condList=[4, 6, 7, 8]
condList

display('extracting cell data')
cd(parentDir)

TextraMeasurements=table('Size',[2,4],'VariableNames',{'columnName','area','valueToCompute','channel'},'VariableTypes',{'string','string','string','string'});
TextraMeasurements.columnName=          {  'avgNuclearDAPI';  'avgCytoplasmDAPI'};
TextraMeasurements.area=                {  'nuclear'       ;  'cytoplasmic'      };
TextraMeasurements.valueToCompute=      {  'average'       ;  'average'         };
TextraMeasurements.channel=             {  'dapi'          ;  'dapi'            };

for i=1:length(condList)
    datetime
    condID=condList(i)
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    %Tcell=extractCells();
    Tcell=extractCells(TextraMeasurements);
 
    % save Tcell in extracted directory
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    
    writetable(Tcell,fullfile(extractedDataDir,'Tcell1.csv')); %for now, save it locally
   
end


%% Extract spot data --> Tspots
% NOTE: maybe should add the gaussian fitting before this, then also extract those
%condList=[1,2,3, 5]
%condList=[5]
%condList=[4, 6, 7, 8]
condList=Tcond.condID
condList
display('------------------------------extracting spots------------------------------')

cd(parentDir)

channelsToExtract=      {'gfp','gfpa','gfpb','gfpc','gfpd','gfpe','tmr','tmra','tmrb','tmrc','tmrd','tmre','alexa','alexaa','alexab','alexac','alexad','cy','cya','cyb','cyc','cyd'};
intensityCutoffs=        [500,  500,    500,   500,   500,   500,   500,  500,    500,   500,   500,   500,   150,    150,     150,     150,     150,     400, 400, 400,  400,  400];
% 20s/image with all channels
%channelsToExtract={'gfpb','gfpc','gfpd','gfpe','tmrb','tmrc','tmrd','tmre','alexab','alexac','alexad','cyb','cyc','cyd'};
% 17s/img for these - not worth it

%intensityCutoffs=repmat(100,1,length(channelsToExtract));


for i=1:length(condList)
    datetime
    condID=condList(i)
    rowTcond=find(Tcond.condID==condID);
    cd(Tcond.Tiff_60X_dir{rowTcond})
    
    % get spots
    Tspots=extractSpots(channelsToExtract,'intensitiesAboveGivenCutoff',intensityCutoffs,'getFittedData',false);
    cd(parentDir)
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    if ~isfolder(extractedDataDir)
        mkdir(extractedDataDir)
    end
    writetable(Tspots,fullfile(extractedDataDir,'Tspots1.csv')); %for now, save it locally
    
end
%% filter out very small and very large cells
% Tcell1 --> Tcell2
% Tspots1 --> Tspots2

%condList=Tcond.condID;

minCellArea=2500;
maxCellArea=50000;

cd(parentDir)
for i=1:length(condList)
    datetime
    condID=condList(i)
    rowTcond=find(Tcond.condID==condID);
    extractedDataDir=Tcond.Extracted_60X_dir{rowTcond};
    
    cd(extractedDataDir)
    Tcell1=readtable('Tcell1.csv');
    Tspots1=readtable('Tspots1.csv');

    % keep only those within
    Tcell2=Tcell1(all([Tcell1.areaCell>=minCellArea,Tcell1.areaCell<=maxCellArea],2),:);    
    Tcell2=mergevars(Tcell2,{'arrayNum','objNum'},'NewVariableName','arrayObjNums');
    
    Tspots1=mergevars(Tspots1,{'arrayNum','objNum'},'NewVariableName','arrayObjNums');
    
    Tspots2=innerjoin(Tspots1,Tcell2(:,{'arrayObjNums'}));
    Tspots2=splitvars(Tspots2,'arrayObjNums','NewVariableNames',{'arrayNum','objNum'});
    
    writetable(Tcell2,'Tcell2.csv')
    writetable(Tspots2,'Tspots2.csv')
    
    
    cd(parentDir)
end

%% select only top N spots per cell for each channel Tspots2 --> Tspots3
% also make Tavg, Tmin, and Tmax, each with all channels in them. 
condList=Tcond.condID
%condList=[1, 2, 3, 5]

channelsToExtract={'gfp','gfpa','gfpb','gfpc','gfpd','gfpe','tmr','tmra','tmrb','tmrc','tmrd','tmre','alexa','alexaa','alexab','alexac','alexad','cy','cya','cyb','cyc','cyd'};
channelBaseNames={'gfp','tmr','alexa','cy'};
spotsPerCell=    [120,  1    , 20     ,80];
channel2numSpotsPerCellMap=containers.Map(channelBaseNames,spotsPerCell)

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/';
cd(parentDir)

TmedIntensity=table('Size',[length(condList),7+length(channelsToExtract)],'VariableTypes',[{'char'},repmat({'double'},1,6+length(channelsToExtract))],'VariableNames',[{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'},channelsToExtract]);
TavgIntensity=table('Size',[length(condList),7+length(channelsToExtract)],'VariableTypes',[{'char'},repmat({'double'},1,6+length(channelsToExtract))],'VariableNames',[{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'},channelsToExtract]);
TminIntensity=table('Size',[length(condList),7+length(channelsToExtract)],'VariableTypes',[{'char'},repmat({'double'},1,6+length(channelsToExtract))],'VariableNames',[{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'},channelsToExtract]);
TmaxIntensity=table('Size',[length(condList),7+length(channelsToExtract)],'VariableTypes',[{'char'},repmat({'double'},1,6+length(channelsToExtract))],'VariableNames',[{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'},channelsToExtract]);


TmedIntensity{:,2:end}=nan;
TavgIntensity{:,2:end}=nan;
TminIntensity{:,2:end}=nan;
TmaxIntensity{:,2:end}=nan;

for i=1:length(condList)
    datetime
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    condName=Tcond.condName{rowTcond};
    cd(Tcond.Extracted_60X_dir{rowTcond})
    
    Tspots2=readtable('Tspots2.csv');
    Tcell2=readtable('Tcell2.csv');
    warning('add back in Tspots1 and Tcell1 read')
    
    numCellThisCond=height(Tcell2);
    

    
    uniqueChannels=unique(Tspots2.channel);
    numUniqueChannels=length(uniqueChannels);
    
    % fill in condition-level info for TavgIntensity and TlowestIntensity
    TmedIntensity(i,{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'})={condName condID numCellThisCond channel2numSpotsPerCellMap('gfp') channel2numSpotsPerCellMap('tmr') channel2numSpotsPerCellMap('alexa') channel2numSpotsPerCellMap('cy')};
    TavgIntensity(i,{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'})={condName condID numCellThisCond channel2numSpotsPerCellMap('gfp') channel2numSpotsPerCellMap('tmr') channel2numSpotsPerCellMap('alexa') channel2numSpotsPerCellMap('cy')};
    TminIntensity(i,{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'})={condName condID numCellThisCond channel2numSpotsPerCellMap('gfp') channel2numSpotsPerCellMap('tmr') channel2numSpotsPerCellMap('alexa') channel2numSpotsPerCellMap('cy')};
    TmaxIntensity(i,{'condName','condID','numCells','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'})={condName condID numCellThisCond channel2numSpotsPerCellMap('gfp') channel2numSpotsPerCellMap('tmr') channel2numSpotsPerCellMap('alexa') channel2numSpotsPerCellMap('cy')};
    
    % go through each channel and fill in Tspots2, TavgIntensity, and TlowestIntensity
    Tspots3=table();
    for iChannel=1:numUniqueChannels
        channel=uniqueChannels{iChannel};
        Tspots2_thisChannel=Tspots2(strcmp(Tspots2.channel,channel),:);
        Tspots2_thisChannel=sortrows(Tspots2_thisChannel,{'intensities'},'descend');
        
        
        % calculate nSpotsToInclude=nCellThisCond*numSpotsPerCell
        indexChannelBaseName=find(ismember(channelBaseNames,{channel,channel(1:end-1)}));
        assert(length(indexChannelBaseName)==1)
        channelBaseName=channelBaseNames{indexChannelBaseName};
        numSpotsPerCell=channel2numSpotsPerCellMap(channelBaseName);
        nSpotsToInclude=round(numCellThisCond*numSpotsPerCell);
        if nSpotsToInclude<200
            warning('low spot count')
        end
        % get the highest-intensity nSpotsToInclude spots, if there are
        % that many spots
        numSpotsTotalThisChannel=height(Tspots2_thisChannel);
        if numSpotsTotalThisChannel>=nSpotsToInclude
            fprintf('condID=%i (%s) channel=%s: numSpotsTotalThisChannel=%i and nSpotsToInclude=%i: GET SPOTS\n',condID,condName,channel,numSpotsTotalThisChannel,nSpotsToInclude)
            % then we have enough spots
            Tspots3_thisChannel=Tspots2_thisChannel(1:nSpotsToInclude,:); %Tspots2 only has top spots
            
            % fill in TavgIntensity and TlowestIntensity and ThighestIntensity
            TmedIntensity{i,{channel}}=median(Tspots3_thisChannel.intensities);
            TavgIntensity{i,{channel}}=mean(Tspots3_thisChannel.intensities);
            TminIntensity{i,{channel}}=min(Tspots3_thisChannel.intensities);
            TmaxIntensity{i,{channel}}=max(Tspots3_thisChannel.intensities);
        else
            Tspots3_thisChannel=Tspots2_thisChannel;
            fprintf('condID=%i (%s) channel=%s: numSpotsTotalThisChannel=%i and nSpotsToInclude=%i: SKIP\n',condID,condName,channel,numSpotsTotalThisChannel,nSpotsToInclude)
            % pad with nan
            nonChannelVariableNames=Tspots3_thisChannel.Properties.VariableNames(~ismember(Tspots3_thisChannel.Properties.VariableNames,'channel'));
            Tspots3_thisChannel(numSpotsTotalThisChannel+1:nSpotsToInclude,nonChannelVariableNames)=array2table(nan(nSpotsToInclude-numSpotsTotalThisChannel,length(nonChannelVariableNames)));
            Tspots3_thisChannel(numSpotsTotalThisChannel+1:nSpotsToInclude,{'channel'})=repmat({''},nSpotsToInclude-numSpotsTotalThisChannel,1);
            % do not fill in TmedIntensity, TavgIntensity, TminIntensity, or TmaxIntensity (leave it nan)
        end
        Tspots3=[Tspots3;Tspots3_thisChannel];
    end
    writetable(Tspots3,'Tspots3.csv')
    writetable(TmedIntensity(i,:),'TmedIntensity.csv')
    writetable(TavgIntensity(i,:),'TavgIntensity.csv')
    writetable(TminIntensity(i,:),'TminIntensity.csv')
    writetable(TmaxIntensity(i,:),'TmaxIntensity.csv')
    
    cd(parentDir)
end

TmedIntensity
TavgIntensity
TminIntensity
TmaxIntensity


%% put tables from all conditions together into one:
% Tcell2, TmedIntensity, TavgIntensity, TminIntensity, TmaxIntensity
% Tspots3 can be put into one but it takes a very long time

condList=Tcond.condID

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/'
extractedDataDirForAllCond='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/extractedData/E159part1rep2_ampclick/';

cd(parentDir)
if ~isfolder(extractedDataDirForAllCond)
        mkdir(extractedDataDirForAllCond)
end
TcellAll=table();
TspotsAll=table();
TmedIntensityAll=table();
TavgIntensityAll=table();
TminIntensityAll=table();
TmaxIntensityAll=table();

for i=1:length(condList)
    condID=condList(i)
    extractedDataDir=Tcond.Extracted_60X_dir{condID};
    cd(extractedDataDir)
    
    Tcell=readtable('Tcell2.csv','ReadVariableNames', true);
    Tspots=readtable('Tspots3.csv','ReadVariableNames', true);
    TmedIntensity=readtable('TmedIntensity.csv');
    TavgIntensity=readtable('TavgIntensity.csv');
    TminIntensity=readtable('TminIntensity.csv');
    TmaxIntensity=readtable('TmaxIntensity.csv');
    
    TcellAll=[TcellAll;[array2table(repmat(condID,height(Tcell),1),'VariableNames',{'condID'}),Tcell]];
    TspotsAll=[TspotsAll;[array2table(repmat(condID,height(Tspots),1),'VariableNames',{'condID'}),Tspots]];
    TmedIntensityAll=[TmedIntensityAll;TmedIntensity];
    TavgIntensityAll=[TavgIntensityAll;TavgIntensity];
    TminIntensityAll=[TminIntensityAll;TminIntensity];
    TmaxIntensityAll=[TmaxIntensityAll;TmaxIntensity];
    
    cd(parentDir)
end

writetable(TcellAll,fullfile(extractedDataDirForAllCond,'Tcell2All.csv'))
writetable(TspotsAll,fullfile(extractedDataDirForAllCond,'Tspots3All.csv'))
writetable(TmedIntensityAll,fullfile(extractedDataDirForAllCond,'TmedIntensityAll.csv'))
writetable(TavgIntensityAll,fullfile(extractedDataDirForAllCond,'TavgIntensityAll.csv'))
writetable(TminIntensityAll,fullfile(extractedDataDirForAllCond,'TminIntensityAll.csv'))
writetable(TmaxIntensityAll,fullfile(extractedDataDirForAllCond,'TmaxIntensityAll.csv'))
TmedIntensityAll
TavgIntensityAll
TminIntensityAll
TmaxIntensityAll



%% create Tcomb - combined table to summarize the channel intensities (median, mean, min, max) but for only 1 'short' and 1 'long' exposure time AND with a normalized median intensity
% also create Tspots4AllNormalized with normalized intensities, based on long/short median scaling at round 6 where short exposures are used
cd(parentDir)

createTspotsAllNormalized=true;

if createTspotsAllNormalized
    TspotsAll=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tspots3All.csv'),'ReadVariableNames', true);
    TspotsAllNormalized=table();
else
    TspotsAll=[];
    TspotsAllNormalized=[];
end

TmedIntensityAll=readtable(fullfile(extractedDataDirForAllCond,'TmedIntensityAll.csv'));
TavgIntensityAll=readtable(fullfile(extractedDataDirForAllCond,'TavgIntensityAll.csv'));
TminIntensityAll=readtable(fullfile(extractedDataDirForAllCond,'TminIntensityAll.csv'));
TmaxIntensityAll=readtable(fullfile(extractedDataDirForAllCond,'TmaxIntensityAll.csv'));
numCond=height(TavgIntensityAll);


Tcomb=table('Size',[numCond,4*15],'VariableTypes',              {'double',    'double',        'char',                  'double',       'double',        'double',        'double',      'char', 'double',                          'double',     'double',     'double',     'double',     'char',        'double',...
                                                                 'double',    'double',        'char',                  'double',       'double',        'double',        'double',      'char', 'double',                          'double',     'double',     'double',     'double',     'char',        'double',...
                                                                 'double',    'double',        'char',                  'double',       'double',        'double',        'double',      'char', 'double',                          'double',     'double',     'double',     'double',     'char',        'double',...
                                                                 'double',    'double',        'char',                  'double',       'double',        'double',        'double',      'char', 'double',                          'double',     'double',     'double',     'double',     'char',        'double'},...
                                                'VariableNames',{'gfp_normMed','gfp_normFactor','gfp_normExpType',      'gfp_shortMed','gfp_shortAvg','gfp_shortMin','gfp_shortMax','gfp_shortName', 'gfp_shortExp',               'gfp_longMed','gfp_longAvg','gfp_longMin','gfp_longMax','gfp_longName','gfp_longExp',...
                                                                 'tmr_normMed','tmr_normFactor','tmr_normExpType',      'tmr_shortMed','tmr_shortAvg','tmr_shortMin','tmr_shortMax','tmr_shortName', 'tmr_shortExp',               'tmr_longMed','tmr_longAvg','tmr_longMin','tmr_longMax','tmr_longName','tmr_longExp',...
                                                                 'alexa_normMed','alexa_normFactor','alexa_normExpType','alexa_shortMed','alexa_shortAvg','alexa_shortMin','alexa_shortMax','alexa_shortName','alexa_shortExp',      'alexa_longMed','alexa_longAvg','alexa_longMin','alexa_longMax','alexa_longName','alexa_longExp',...
                                                                 'cy_normMed','cy_normFactor','cy_normExpType',         'cy_shortMed','cy_shortAvg','cy_shortMin','cy_shortMax','cy_shortName',        'cy_shortExp',             'cy_longMed','cy_longAvg','cy_longMin','cy_longMax','cy_longName','cy_longExp'});

% set default to nan
for iCol=1:width(Tcomb)
    colName=Tcomb.Properties.VariableNames{iCol};
    colVect=Tcomb.(colName);
    if isnumeric(colVect)
        Tcomb.(colName)=nan(height(Tcomb),1);
    end
end
                                            
Ttemp=TmedIntensityAll(:,{'condName','condID','numCells'});
Ttemp=join(Ttemp,Tcond(:,{'condName','AmpRound','SpecialClick'}),'Keys','condName');
Ttemp=join(Ttemp,TmedIntensityAll(:,{'condName','gfp_spotsPerCell','tmr_spotsPerCell','alexa_spotsPerCell','cy_spotsPerCell'}),'Keys','condName');
Tcomb=[Ttemp,Tcomb];
Ttemp=[];

ampRoundToNormalizeTo=1;
firstAmpRoundToUseShortExp=6; % not applicable for No Click conditions

rowOfAmpRoundToNormalizeTo=find(Tcomb.AmpRound==ampRoundToNormalizeTo);
assert(isscalar(rowOfAmpRoundToNormalizeTo))
rowOfFirstAmpRoundToUseShortExp=find(Tcomb.AmpRound==firstAmpRoundToUseShortExp);
assert(length(rowOfFirstAmpRoundToUseShortExp)==1)
%TavgCombinedExposures=[array2table(cell(numCond,1),'VariableNames',{'condName'},...
%                       array2table(nan(numCond,9),'VariableNames',{'condID','gfp_avg','gfp_exposure','tmr_avg','tmr_exposure','alexa_avg','alexa_exposure','cy_avg','cy_exposure'})];


channelBaseNames={'gfp','tmr','alexa','cy'};
longExp=        [1000,  1000, 500,   500];
shortExp=        [100,  100, 100,   100];

TchannelExpToUse=[array2table(channelBaseNames','VariableNames',{'channelBaseName'}),...
                array2table(longExp','VariableNames',{'longExp'}),...
                array2table(shortExp','VariableNames',{'shortExp'})];

channelExpToLastLetterMap=containers.Map([ 2, 5,  50, 100,500,1000],...
                                         {'','a','b','c','d','e'} );


for i=1:numCond
    condID=Tcomb.condID(i);
    ampRound=Tcomb.AmpRound(i);
    SpecialClick=Tcomb.SpecialClick{i};
    
    for iChannel=1:height(TchannelExpToUse)
        % what are channel baseName, short exposure, long exposure, short channelname
        channel_baseName=TchannelExpToUse.channelBaseName{iChannel};
        channel_shortExp=TchannelExpToUse.shortExp(strcmp(TchannelExpToUse.channelBaseName,channel_baseName));
        channel_longExp=TchannelExpToUse.longExp(strcmp(TchannelExpToUse.channelBaseName,channel_baseName));
        channel_shortName=[channel_baseName,channelExpToLastLetterMap(channel_shortExp)];
        channel_longName=[channel_baseName,channelExpToLastLetterMap(channel_longExp)];
        
        % get channel_shortMed, channel_shortAvg, channel_shortMin, channel_shortMax
        channel_shortMed=TmedIntensityAll.(channel_shortName)(TmedIntensityAll.condID==condID);
        channel_shortAvg=TavgIntensityAll.(channel_shortName)(TavgIntensityAll.condID==condID);
        channel_shortMin=TminIntensityAll.(channel_shortName)(TminIntensityAll.condID==condID);
        channel_shortMax=TmaxIntensityAll.(channel_shortName)(TmaxIntensityAll.condID==condID);
        
        % get channel_longAvg, channel_longMin, channel_longMax
        channel_longMed=TmedIntensityAll.(channel_longName)(TmedIntensityAll.condID==condID);
        channel_longAvg=TavgIntensityAll.(channel_longName)(TavgIntensityAll.condID==condID);
        channel_longMin=TminIntensityAll.(channel_longName)(TminIntensityAll.condID==condID);
        channel_longMax=TmaxIntensityAll.(channel_longName)(TmaxIntensityAll.condID==condID);
        
        
        % put the above in TavgComb
        Tcomb.([channel_baseName,'_shortExp'])(i)=channel_shortExp;
        Tcomb.([channel_baseName,'_longExp'])(i)=channel_longExp;
        Tcomb.([channel_baseName,'_shortName']){i}=channel_shortName;
        Tcomb.([channel_baseName,'_longName']){i}=channel_longName;
        
        Tcomb.([channel_baseName,'_shortMed'])(i)=channel_shortMed;
        Tcomb.([channel_baseName,'_shortAvg'])(i)=channel_shortAvg;
        Tcomb.([channel_baseName,'_shortMin'])(i)=channel_shortMin;
        Tcomb.([channel_baseName,'_shortMax'])(i)=channel_shortMax;
        
        Tcomb.([channel_baseName,'_longMed'])(i)=channel_longMed;
        Tcomb.([channel_baseName,'_longAvg'])(i)=channel_longAvg;
        Tcomb.([channel_baseName,'_longMin'])(i)=channel_longMin;
        Tcomb.([channel_baseName,'_longMax'])(i)=channel_longMax;
       
        % for this channel, what is the first round's long-exposure time median intensity for normalization
        channel_longMedOfRoundToNormalizeTo=TmedIntensityAll.(channel_longName)(rowOfAmpRoundToNormalizeTo);
        
        % fill in the channel_normMed
        if or(ampRound<firstAmpRoundToUseShortExp,~isempty(SpecialClick)) % use long exposure
            channel_normFactorExpTime=1;
            channel_normFactor=channel_normFactorExpTime * (1/channel_longMedOfRoundToNormalizeTo);
            
            channel_normMed=channel_normFactor * channel_longMed;
            normExpType='long';
            
            channelNameOfNormalizedData=channel_longName;
            exposureForNormalizedData=channel_longExp;
        else % use short exposure
            % get norm factor from firstAmpRoundToUseShortExp
            channel_normFactorExpTime=Tcomb.([channel_baseName,'_longMed'])(rowOfFirstAmpRoundToUseShortExp)/Tcomb.([channel_baseName,'_shortMed'])(rowOfFirstAmpRoundToUseShortExp);
            channel_normFactor=channel_normFactorExpTime * (1/channel_longMedOfRoundToNormalizeTo);
            
            % calculate normalized intensity
            channel_normMed=channel_normFactor * channel_shortMed;
            normExpType='short';
            
            channelNameOfNormalizedData=channel_shortName;
            exposureForNormalizedData=channel_shortExp;
        end

            
        Tcomb.([channel_baseName,'_normFactor'])(i)=channel_normFactor;
        Tcomb.([channel_baseName,'_normMed'])(i)=channel_normMed;
        Tcomb.([channel_baseName,'_normExpType']){i}=normExpType;
        
        % now create TspotsAllNormalized
        if createTspotsAllNormalized
            TspotsThisCondThisChannelNorm=TspotsAll(all([TspotsAll.condID==condID,strcmp(TspotsAll.channel,channelNameOfNormalizedData)],2),:);
            nSpotsTCC=height(TspotsThisCondThisChannelNorm);
            
            TspotsThisCondThisChannelNorm.channelBaseName=repmat({channel_baseName},nSpotsTCC,1);
            TspotsThisCondThisChannelNorm.expTime=repmat(exposureForNormalizedData,nSpotsTCC,1);
            TspotsThisCondThisChannelNorm.normFactor=repmat(channel_normFactor,nSpotsTCC,1);
            TspotsThisCondThisChannelNorm.intensitiesNorm= channel_normFactor * TspotsThisCondThisChannelNorm.intensities;
        
            TspotsAllNormalized=[TspotsAllNormalized;TspotsThisCondThisChannelNorm];
        end
        
        
    end
    
end

% add in AmpFactorFromPrevious
rowsForAmpFactorFromPrevious=[1 2 3 4 5 6];
normDataVariables=Tcomb.Properties.VariableNames(endsWith(Tcomb.Properties.VariableNames,'_normAvg'));
for iCol=1:length(normDataVariables)
    normDataVariable=normDataVariables{iCol};
    channelBaseName=normDataVariable(1:regexp(normDataVariable,'_')-1);
    foldChangeFromPrevious=nan(height(Tcomb),1);
    foldChangeFromPrevious(rowsForAmpFactorFromPrevious(2:end))=(Tcomb.(normDataVariable)(rowsForAmpFactorFromPrevious(2:end)))./(Tcomb.(normDataVariable)(rowsForAmpFactorFromPrevious(1:end-1)));
    Tcomb=addvars(Tcomb,foldChangeFromPrevious,'After',normDataVariable,'NewVariableNames', [channelBaseName,'_foldChangeFromPrevious']);
end

cd(parentDir)
% create TavgCombConcise
%conciseVariableNames=TavgComb.Properties.VariableNames(~endsWith(TavgComb.Properties.VariableNames,{'shortAvg','shortName','longAvg','longName'}));
conciseVariableNames=Tcomb.Properties.VariableNames(~endsWith(Tcomb.Properties.VariableNames,{'shortAvg','shortName','longAvg','longName','normFactor','normExpType','shortExp','longExp'}));

TcombConcise=Tcomb(:,conciseVariableNames);

Tcomb
TcombConcise

writetable(Tcomb,fullfile(extractedDataDirForAllCond,filesep,'Tcomb.csv'));
writetable(TcombConcise,fullfile(extractedDataDirForAllCond,filesep,'TcombConsise.csv'));

% add on to TspotsAllNormalized: condName     AmpRound    SpecialClick
TspotsAllNormalized=join(TspotsAllNormalized,Tcond(:,{'condID','condName','AmpRound','SpecialClick'}),'Keys',{'condID'});
TspotsAllNormalized=movevars(TspotsAllNormalized,{'condName','AmpRound','SpecialClick'},'After','condID');
writetable(TspotsAllNormalized,fullfile(extractedDataDirForAllCond,filesep,'Tspots4AllNormalized.csv'));
%% inspect thresholds by setting them in the actual data files then using thresholdGUI to inspect
% cd(parentDir)
% %extractedDataDirForAllCond='extractedData/E159_scan_amp/analysis1_amp';
% extractedDataDirForAllCond='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper/extractedData/E159part1rep2_ampclick/';
% Tcomb=readtable(fullfile(extractedDataDirForAllCond,filesep,'Tcomb.csv'))
% %% View thresholds in thresholdGUI
% % condList=[3];
% % channels={'gfpd','tmrd','alexad','cyc'};
% % thresholds=[2048.4 1678.4 674.26 387.83]; %s
% 
% condList=[1];
% channels={'gfpe','tmre','alexad','cyd'}; %
% thresholds=[610.5 558.95 269.06 707.99];%
% 
% %arrayNum=1;
% %objNum=1;
% parentDir='/Volumes/IAND_08';
% cd(parentDir)
% for i=1:length(condList)
%     datetime
%     condID=condList(i)
%     rowTcond=find(Tcond.condID==condID);
%     cd(Tcond.Tiff_60X_dir{rowTcond})
%     
% 	%tools=improc2.launchImageObjectTools;
% 	%tools.navigator.tryToGoToArray(arrayNum)
% 	%tools.navigator.tryToGoToObj(objNum)
%     %setThreshold(channels,thresholds,tools) % set just the object selected in tools
%     
%     setThreshold(channels,thresholds) % set all objects in current folder
%     
%     
%     
%     improc2.launchThresholdGUI
%     
%     cd(parentDir)
% end
% %%
% %tools=improc2.launchThresholdGUI;
% %tools.browsingTools.navigator.tryToGoToArray(1)
% %tools.browsingTools.navigator.tryToGoToObj(1)
% %tools.rnaChannelSwitch.setChannelName('tmre')
