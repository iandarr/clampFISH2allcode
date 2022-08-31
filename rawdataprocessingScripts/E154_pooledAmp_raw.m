% E154_pooledAmp_raw
% pooled amplification
%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E154_pooled_amp_Tcond.xlsx');


%% Nd2-->Tiff, 60X data (p2 1to125)
inDirFirstRun='data/E154_pooled_amp/Nd2/s6_60X_p2_1to125';
inDirSecondRun='data/E154_pooled_amp/Nd2/s6.2 60X p2 moreData';
listND2_inDirFirstRun=listFilesWithExtension('.nd2','directoryPath',inDirFirstRun);
listND2_inDirSecondRun=listFilesWithExtension('.nd2','directoryPath',inDirSecondRun);

indices=find(strcmp(Tcond.amplifierDilution,'1to125'))';
%indices=[38 43]; %
for i=indices
    condName=Tcond.condName{i};
    outDir=Tcond.Tiff_60X_dir{i};
    mkdir(outDir)
    
    % get both (or one) ND2 file and export to Tiff
    %get firstRun file path
    filenametext_FirstRun=Tcond.Nd2_60X_filenametext_FirstRun(i);
    listInd_FirstRun=find(contains(listND2_inDirFirstRun,filenametext_FirstRun));
    assert(length(listInd_FirstRun)==1)
    filepath_FirstRun=[inDirFirstRun,filesep,listND2_inDirFirstRun{listInd_FirstRun}];
    
    %get secondRun file path, if it exists
    filenametext_SecondRun=Tcond.Nd2_60X_filenametext_SecondRun(i);
    listInd_SecondRun=find(contains(listND2_inDirSecondRun,filenametext_SecondRun));
    if ~isnan(listInd_SecondRun)
        assert(length(listInd_SecondRun)==1)
        filepath_SecondRun=[inDirSecondRun,filesep,listND2_inDirSecondRun{listInd_SecondRun}];
        filepath_AllRuns={filepath_FirstRun; filepath_SecondRun};
    else
        filepath_SecondRun=nan;
        filepath_AllRuns={filepath_FirstRun};
    end
    
    % output Tiffs from firstRun
    %nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
    % output Tiffs from firstRun and secondRun (if there is one)
    nd2toTiff(filepath_AllRuns,'outDir',outDir,'outputXYpositions',true)
    
end
cd(parentDir)
%% Nd2-->Tiff, 60X data (p1 1to250)
% inDirFirstRun='data/E154_pooled_amp/Nd2/s5_60X_p1_1to250';


%% Nd2-->Tiff, 20X data (p1 1to125)
inDirFirstRun='data/E154_pooled_amp/Nd2/s2 20X p1';
listND2_inDirFirstRun=listFilesWithExtension('.nd2','directoryPath',inDirFirstRun);

indices=find(strcmp(Tcond.amplifierDilution,'1to250'))';

for i=indices
    condName=Tcond.condName{i};
    outDir=Tcond.Tiff_20X_dir{i};
    mkdir(outDir)
    
    % get both (or one) ND2 file and export to Tiff
    %get firstRun file path
    filenametext_FirstRun=Tcond.Nd2_20X_filenametext(i);
    listInd_FirstRun=find(contains(listND2_inDirFirstRun,filenametext_FirstRun));
    assert(length(listInd_FirstRun)==1)
    filepath_FirstRun=[inDirFirstRun,filesep,listND2_inDirFirstRun{listInd_FirstRun}];
    
    % output Tiffs from firstRun
    %nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
    % output Tiffs from firstRun and secondRun (if there is one)
    nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
end
cd(parentDir)
%% Nd2-->Tiff, 20X data (p2 1to125)
inDirFirstRun='data/E154_pooled_amp/Nd2/s4 20X p2';
listND2_inDirFirstRun=listFilesWithExtension('.nd2','directoryPath',inDirFirstRun);

indices=find(strcmp(Tcond.amplifierDilution,'1to125'))';
condList=Tcond.condID(strcmp(Tcond.amplifierDilution,'1to125'))';
for i=1:length(condList)
    condID=condList(i);
    rowTcond=find(Tcond.condID==condID);
    
    condName=Tcond.condName{rowTcond};
    outDir=Tcond.Tiff_20X_dir{rowTcond};
    mkdir(outDir)
    
    % get both (or one) ND2 file and export to Tiff
    %get firstRun file path
    filenametext_FirstRun=Tcond.Nd2_20X_filenametext(i);
    listInd_FirstRun=find(contains(listND2_inDirFirstRun,filenametext_FirstRun));
    assert(length(listInd_FirstRun)==1)
    filepath_FirstRun=[inDirFirstRun,filesep,listND2_inDirFirstRun{listInd_FirstRun}];
    
    % output Tiffs from firstRun
    %nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
    % output Tiffs from firstRun and secondRun (if there is one)
    nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
end
cd(parentDir)