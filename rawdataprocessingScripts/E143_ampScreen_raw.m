% E143_ampScreen_raw
% amplifier screen

%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E143_amplifierScreen_Tcond.xlsx');

%% Nd2-->Tiff, 60X data
outParentDir=['Tiff_data',filesep,'60X',filesep'];
inDirFirstRun=['Nd2_raw_data',filesep,'60X_FirstRun'];
listND2_inDirFirstRun=listFilesWithExtension('.nd2','directoryPath',inDirFirstRun);
listND2_inDirSecondRun=listFilesWithExtension('.nd2','directoryPath',inDirSecondRun);

for i=1:height(Tcond)
    condName=Tcond.condName{i};
    outDir=Tcond.Tiff_60X_dir{i}
    mkdir(outDir)
    
    % get both (or one) ND2 file and export to Tiff
    % note condID=42 (EGFR_ser12) only has 60X first run
    
    %get firstRun file path
    filenametext_FirstRun=Tcond.Nd2_60X_filenametext_FirstRun(i);
    listInd_FirstRun=find(contains(listND2_inDirFirstRun,filenametext_FirstRun));
    assert(length(listInd_FirstRun)==1)
    filepath_FirstRun=[inDirFirstRun,filesep,listND2_inDirFirstRun{listInd_FirstRun}];
    
    %     %get secondRun file path
    %     filenametext_SecondRun=Tcond.Nd2_60X_filenametext_SecondRun(i);
    %     listInd_SecondRun=find(contains(listND2_inDirSecondRun,filenametext_SecondRun));
    %     assert(length(listInd_SecondRun)==1)
    %     filepath_SecondRun=[inDirSecondRun,filesep,listND2_inDirSecondRun{listInd_SecondRun}];
    
    nd2toTiff(filepath_FirstRun,'outDir',outDir,'outputXYpositions',true)
    
end

%% Nd2-->Tiff, 20X data
inDirFirstRun=['Nd2_raw_data',filesep,'20X'];
listND2_inDir=listFilesWithExtension('.nd2','directoryPath',inDirFirstRun);

listOutDir=Tcond.Tiff_20X_dir;
listFilenametext=Tcond.Nd2_20X_filenametext;
%condIDAlreadyExported=[1 2 3 4 11 15 16 26 30 31 41 45 47 48]';
%condIDNeedToExport=Tcond.condID(~ismember(Tcond.condID,condIDAlreadyExported))';

for i=1:height(Tcond)
    condName=Tcond.condName{i};
    outDir=listOutDir{i};
    mkdir(outDir)

    % get both (or one) ND2 file and export to Tiff
    % note condID=42 (EGFR_ser12) only has 60X first run

    %get firstRun file path
    filenametext=listFilenametext(i);
    listInd=find(contains(listND2_inDir,filenametext));
    assert(length(listInd)==1)
    filepath=[inDirFirstRun,filesep,listND2_inDir{listInd}];

    nd2toTiff(filepath,'outDir',outDir,'outputXYpositions',true)

end