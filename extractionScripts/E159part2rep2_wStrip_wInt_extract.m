% E159part2rep2_scanWithStrip_withMeanInt_extract
% first run E159part2rep2_tilebg before running this
%   scanID=3 (part of well A1, drug-naive cells,     with extra imaging rounds (R13rep, R1s, R2s, R3s, R4) with post-strip data (R1s,R2s,R3s) and repeat of R1 data (R4)),

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tscan.xlsx'));
extractedDataDirIn=fullfile(parentDir,'paper/extractedData/E159part2rep2_scan/scan3_withStripAndRepeatData',filesep);
extractedDataDirOut=fullfile(parentDir,'paper/extractedData/E159part2rep2_scan/scan3_withStripAndRepeatData_withMeans',filesep);
cd(parentDir)
%scanIDsList=[7,8]';
scanIDsList=3;
Tchannels=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scanWithStrip_ChannelDesc.xlsx'));
Tchannels=Tchannels(~ismember(Tchannels.channelName,{'Merged_Rcyto_YFP_DAPI','Merged_Rcyto_YFP_DAPI_cp_masks','R123_maxYFP'}),:)
%preStitchedScanFilelist=join([Tchannels.channelName,repmat({'.tif'},height(Tchannels),1)],'',2);

backgroundTilesDir='/Volumes/IAND_08/rawData/E159part2rep2_scan/scan/scan3_WellA3withStrip/BackgroundTiles';
numChannels=height(Tchannels);

imagingCycleList={'R1','R1rep','R1s','R2','R2s','R3','R3s','R4','Rcyto'};
numImagingCycles=length(imagingCycleList);
Tchannels.imagingCycle=extractBefore(Tchannels.channelName,'_');

%% Make TcellMeanIntensity.csv: with mean fluorescence intensity
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

    %warning('changed this to start at subregion 14')
    for iSubregion=1:numSubregion
        subregionDir=subregionDirList{iSubregion};
        fprintf('scanID=%i: getting mean intensities in cells for subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))

        % Get the average intensity of each channel
        maskPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,'Merged_Rcyto_YFP_DAPI_cp_masks.tif');
        maskImg=imread(maskPath);

        % insert getting average intensity for all channels here!
        numCells=max(maskImg,[],'all');
        TcellMeanIntens=array2table(nan(numCells,numChannels),'VariableNames',Tchannels.channelName');
        TcellMeanIntens.cellID=[1:numCells]';
        TcellMeanIntens=movevars(TcellMeanIntens,'cellID','Before',1);
        %channelList=
        fprintf('  iChannel working on= ')
        for iChannel=1:numChannels
            fprintf('%i ',iChannel)
            channelName=Tchannels.channelName{iChannel};
            channelImgFilename=[channelName,'.tif'];
            channelImgPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,channelImgFilename);
            channelImg=imread(channelImgPath);

            % mean intensities of image
            TpropsChannel=regionprops("table",maskImg,channelImg,{'MeanIntensity'});
            TcellMeanIntens.(channelName)=TpropsChannel{:,{'MeanIntensity'}};

        end
        fprintf('\n')
        fprintf('   writing table to file\n')
        outTablePath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanIntensity.csv');
        writetable(TcellMeanIntens,outTablePath)


    end
end

%% Make TcellMeanBackground.csv: with mean fluorescence intensity of background
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
        fprintf('scanID=%i for subregionDir=%s: getting mean background intensities in cell segmentation areas for starting %s\n',scanID,subregionDir,char(datetime))

        % Get the average background intensity
        maskPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,'Merged_Rcyto_YFP_DAPI_cp_masks.tif');
        maskImg=imread(maskPath);

        % insert getting average intensity for all channels here!
        numCells=max(maskImg,[],'all');
        TcellMeanBackground=array2table(nan(numCells,numChannels),'VariableNames',Tchannels.channelName');
        TcellMeanBackground.cellID=[1:numCells]';
        TcellMeanBackground=movevars(TcellMeanBackground,'cellID','Before',1);


        for iCycle=1:numImagingCycles
            imagingCycle=imagingCycleList{iCycle};
            channelNamesThisCycle=Tchannels.channelName(strcmp(Tchannels.imagingCycle,imagingCycle));
            numChannelsThisCycle=length(channelNamesThisCycle);
            % get image coordinates of the cell in the original tiles
            rowsImgPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,[imagingCycle,'_ImgRow.tif']);
            rowsImg=imread(rowsImgPath);
            Trows=regionprops("table",maskImg,rowsImg,{'PixelValues'});
            colsImgPath=fullfile(subregionParentDir,filesep,subregionDir,filesep,[imagingCycle,'_ImgCol.tif']);
            colsImg=imread(colsImgPath);
            Tcols=regionprops("table",maskImg,colsImg,{'PixelValues'});

            fprintf('  imagingCycle=%8s  iChannel working on= ',imagingCycle)
            for iChannel=1:numChannelsThisCycle
                fprintf('%i ',iChannel)
                channelName=channelNamesThisCycle{iChannel};
                backgroundTileFilename=[extractBefore(channelName,'ms'),'ms.tif'];
                backgroundTilePath=fullfile(backgroundTilesDir,filesep,backgroundTileFilename);
                backgroundTile=imread(backgroundTilePath);
                sizeBgTile=size(backgroundTile);
                meanBackgroundIntensities=nan(numCells,1);
                for iCell=1:numCells
                    ind=sub2ind(sizeBgTile,Trows.PixelValues{iCell},Tcols.PixelValues{iCell});
                    meanBackgroundIntensities(iCell)=mean(backgroundTile(ind));
                end
                TcellMeanBackground.(channelName)=meanBackgroundIntensities;

            end
            fprintf('\n')
        end
        fprintf('   writing table to file\n')
        outTablePath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanBgInt.csv');
        writetable(TcellMeanBackground,outTablePath)


    end
end

%% TcellMeanSignalInt = Subtraction: TcellMeanCellInt - TcellMeanBgInt
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

    %warning('changed this to start at subregion 14')
    for iSubregion=1:numSubregion
        subregionDir=subregionDirList{iSubregion};
        fprintf('scanID=%i: calculating TcellMeanSignalInt for subregionDir=%s starting %s\n',scanID,subregionDir,char(datetime))

        % get raw cell intensity and Background intensity tables
        meanCellIntPath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanCellInt.csv');
        rawInt=readtable(meanCellIntPath);
        bgIntPath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanBgInt.csv');
        bgInt=readtable(bgIntPath);
        assert(isequal(rawInt.Properties.VariableNames,bgInt.Properties.VariableNames))
        assert(isequal(rawInt.cellID,bgInt.cellID))

        % calculate signal (raw cell intensity minus background)
        sigInt=array2table([rawInt{:,1},rawInt{:,2:end} - bgInt{:,2:end}],'VariableNames',rawInt.Properties.VariableNames);
        % append raw_ or bg_ or sig_ to VariableNames
        rawInt.Properties.VariableNames(2:end)=join([repmat({'raw_'},width(rawInt)-1,1),[rawInt.Properties.VariableNames(2:end)]'],'',2)';
        bgInt.Properties.VariableNames(2:end)=join([repmat({'bg_'},width(bgInt)-1,1),[bgInt.Properties.VariableNames(2:end)]'],'',2)';
        sigInt.Properties.VariableNames(2:end)=join([repmat({'sig_'},width(sigInt)-1,1),[sigInt.Properties.VariableNames(2:end)]'],'',2)';
        fprintf('   writing table for mean: signal (raw cell - background), raw, and background to file\n')
        % combine into a single table TcellMeanCellSigIntBg
        TcellMeanCellSigRawBgCombined=[sigInt,rawInt(:,2:end),bgInt(:,2:end)];
        outTablePath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanCellSigRawBg.csv');
        writetable(TcellMeanCellSigRawBgCombined,outTablePath)


    end
end

%% combine each subregion's TcellMeanCellSigRawBg into one
scanIDsList=3;

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
        %cd([subregionParentDir,filesep,subregionDir])
        %opts = detectImportOptions('Tcell2.csv');
        %opts = setvartype(opts, {'segID'}, 'uint32');
        %Tcell=readtable('Tcell2.csv',opts);
        TcellMeansPath=fullfile([subregionParentDir,filesep,subregionDir,filesep],'TcellMeanCellSigRawBg.csv');
        Tcell=readtable(TcellMeansPath);
        Tcell.segID=Tcell.cellID;
        Tcell=movevars(Tcell,'segID','Before',1);
        Tcell.cellID=[]; % what I called cellID should really be segID, since cellID is for the whole scan not just the subregion

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
        %cd(parentDir)
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

fullPathName=fullfile(extractedDataDirOut,'TcellAllMeanSigRawBg.csv');
writetable(TcellAll,fullPathName)

%% TcellAll1_withMeans: Merge of TcellAll1 and TcellMeanCellSigRawBg
TcellAll1=readtable(fullfile(extractedDataDirIn,'TcellAll1.csv')); % with spot counts
TcellMeanCellSigRawBg=readtable(fullfile(extractedDataDirOut,'TcellAllMeanSigRawBg.csv')); % with mean intensities
% isResistant is last - move it to 9th position after segID.
TcellAll1=movevars(TcellAll1,'isResistant','After','segID');
TcellMeanCellSigRawBg=movevars(TcellMeanCellSigRawBg,'isResistant','After','segID');
assert(isequal(TcellMeanCellSigRawBg(:,1:9),TcellAll1(:,1:9)))
TcellAll1_WithMeans=[TcellAll1,TcellMeanCellSigRawBg(:,10:end)];
fullPathName=fullfile(extractedDataDirOut,'TcellAll1_WithMeans.csv');
writetable(TcellAll1_WithMeans,fullPathName)

%% Apply QC filters TcellAll1_withMeans.csv --> TcellAll1_withMeans_withQC.csv and TcellAll2_withMeans.csv
% considering these readout cycles ('rounds') for UBC counts: R1, R2, R3

TcellAll=readtable(fullfile(extractedDataDirOut,'TcellAll1_withMeans.csv'));
%TcellAll.isResistant=logical(TcellAll.isResistant);
%TcellAll=TcellAll(TcellAll.isResistant,:);
%warning('remove ~isResistant filter')

% QC critera 1: only subregions with meanUBC > = subregionMinMeanUBC
Tsr=groupsummary(TcellAll(:,{'uniqueSubregionInd','subregionInd','subregionRow','subregionCol','scanID','R1_YFP_100ms_UBC','R2_YFP_100ms_UBC','R3_YFP_100ms_UBC'}),'uniqueSubregionInd','mean');
Tsr=renamevars(Tsr,{'mean_subregionInd','mean_subregionRow','mean_subregionCol','mean_scanID'},{'subregionInd','subregionRow','subregionCol','scanID'});

subregionMinMeanUBC=25; % min of the average of all UBC per cell for each subregion. If you masked a lot of spots (and their area is now zero) then this would produce a lot of zero-count cells, and would bring down this average
% ^ would probably need to be increased for well B3 (scan 6,8), but this was imaged in a smaller region in the middle of the well, so coverslip edge  problem not applicable

Tsr.passSubregion=all([Tsr.mean_R1_YFP_100ms_UBC,Tsr.mean_R2_YFP_100ms_UBC,Tsr.mean_R3_YFP_100ms_UBC]>=subregionMinMeanUBC,2);
TcellAll=join(TcellAll,Tsr(:,{'uniqueSubregionInd','passSubregion'}));

% QC criteria 2: For all 3 rounds, UBC count >= minUBC
%   Note: cells can have zero UBC if the spots were masked in dentist2, or
%   if truly no spots above threshold were in that segmentation

UBC3Rounds=TcellAll{:,{'R1_YFP_100ms_UBC','R2_YFP_100ms_UBC','R3_YFP_100ms_UBC'}};

minUBC=4; % per cell;
passMinUBC=all(UBC3Rounds>=minUBC,2);
TcellAll.passMinUBC=passMinUBC;

% QC criteria 3: UBC >= minUBCPerUm2 * CellArea
%   Note: larger cells have more UBC. If it doesn't, it's likely out of focus
minUBCPerUm2=0.025;

passMinUBCPerArea = min(UBC3Rounds,[],2)>= minUBCPerUm2 * TcellAll.Area_Um2;
TcellAll.passMinUBCPerArea=passMinUBCPerArea;

% QC critera 4: Percent UBC error cutoff
minPcrErrUBC=50; % 50 = UBC from all rounds should have less than or equal to 50% difference from median of the 5 rounds
passPctUBC=all((100 * abs(UBC3Rounds./median(UBC3Rounds,2) - 1)) <= minPcrErrUBC,2);
TcellAll.passPctUBC=passPctUBC;

% passes all QC checks
passAllQC=all([TcellAll.passSubregion, TcellAll.passMinUBC, TcellAll.passMinUBCPerArea, TcellAll.passPctUBC],2);
TcellAll.passAllQC=passAllQC;

numOrig=height(TcellAll);
numPassedAllQC=sum(TcellAll.passAllQC);
numPassedSubregion=sum(TcellAll.passSubregion);
numGoodSubregion=sum(Tsr.passSubregion);
numTotSubregion=height(Tsr);

fprintf(['QC step 1 (subregion-level): with subregionMinMeanUBC=%.2f, %i subregions (of %i) passed (%i of %i cells, or %.2f%%)\n',...
         'QC step 2 (cell-level):      %i of %i cells (%.2f%%) passed cell-specific QC checks \n',...
         'Overall: %i (%.2f%%) of cells passed all QC\n'],...
                           subregionMinMeanUBC, numGoodSubregion,numTotSubregion, numPassedSubregion, numOrig, 100*numPassedSubregion/numOrig,...
                           numPassedAllQC,numPassedSubregion,100*numPassedAllQC/numPassedSubregion,...
                           numPassedAllQC,100*numPassedAllQC/numOrig)

%groupsummary(TcellAll(:,{'subregionInd','R1_YFP_UBC'}),'subregionInd','mean')
%groupsummary(TcellAll(:,{'subregionInd','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC','passMinUBC','passPctUBC','passAllQC'}),'subregionInd','mean')
%groupsummary(TcellAll(:,{'uniqueSubregionInd','scanID','subregionInd','R1_YFP_UBC','R2_YFP_UBC','R3_YFP_UBC','passMinUBC','passPctUBC','passAllQC'}),'uniqueSubregionInd','mean')

writetable(TcellAll,fullfile(extractedDataDirOut,'TcellAll1_withMeans_withQC.csv'))

TcellAll2=TcellAll(TcellAll.passAllQC,:); % just what passes QC
writetable(TcellAll2,fullfile(extractedDataDirOut,'TcellAll2_withMeans.csv'))

