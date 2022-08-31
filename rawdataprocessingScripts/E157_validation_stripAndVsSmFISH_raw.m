% E157_VsSmFISH_raw
% part1: readout probe stripping
% part2: clampFISH 2.0 vs. smFISH correlation
%% Table of all conditions

%parentDir='/Volumes/IAND_04';
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E157_validation_Tcond.xlsx')
if ~all(Tcond.condID==[1:height(Tcond)]')
    error('condID must be equivalent its row number in this code')
end

%% Nd2-->Tiff, part1, 20X data via nd2ToTiff
cd(parentDir)

%indices=find(strcmp(Tcond.purpose,'strip'))';
indices=1:10;

%magnificationsToExportTiffsFor={'60X','20X','10X'};
%magnificationsToExportTiffsFor={'20X','10X'};
magnificationsToExportTiffsFor={'20X'};
for i=1:length(magnificationsToExportTiffsFor)
    magnification=magnificationsToExportTiffsFor{i};

    fprintf('%%%%% magnification = %s %%%%%\n',magnification)


for ii=indices
    datetime
    condName=Tcond.condName{ii}

    % get both (or one) ND2 file and export to Tiff
    %get firstRun file path
    Nd2_dir=Tcond.(['Nd2_',magnification,'_dir']){ii};

    filenametext=Tcond.(['Nd2_',magnification,'_filenametext'])(ii);

    listND2_inDir=listFilesWithExtension('.nd2','directoryPath',Nd2_dir);
    listInd=find(all([contains(listND2_inDir,filenametext),~startsWith(listND2_inDir,'.')],2));

    if length(listInd)==1
    %assert(length(listInd)==1)
    filepath=[Nd2_dir,filesep,listND2_inDir{listInd}];

    if isfile(fullfile(Nd2_dir,'channelMapInfo.xlsx'))
        useCustomChannelMap=true;
        channelMapTable=readtable(fullfile(Nd2_dir,'channelMapInfo.xlsx'));
        channelMap=containers.Map(channelMapTable{:,1}',channelMapTable{:,2}');
    else
        useCustomChannelMap=false;
    end

    % create directory for Tiffs
    outDir=Tcond.(['Tiff_',magnification,'_dir']){ii}
    if isdir(outDir)
        warning('overwriting data in %s',outDir)
        rmdir(outDir,'s')
        mkdir(outDir)
    else
        mkdir(outDir)
    end

    % output Tiffs
    if useCustomChannelMap
        nd2toTiff(filepath,'outDir',outDir,'outputXYpositions',true,'channelMap',channelMap)
    else
        nd2toTiff(filepath,'outDir',outDir,'outputXYpositions',true)
    end

    else
        warning('skipping condName=%s Mag=%s, did not find one file with filenametext in folder %s',condName,magnification,Nd2_dir)
    end
end
end
cd(parentDir)

%% Nd2-->Tiff, part2, 60X,20X,and 10X data via nd2ToTiff
cd(parentDir)

%indices=find(strcmp(Tcond.purpose,'strip'))';
indices=11:22;

magnificationsToExportTiffsFor={'60X','20X','10X'};
%magnificationsToExportTiffsFor={'20X','10X'};
%magnificationsToExportTiffsFor={'20X'};
for i=1:length(magnificationsToExportTiffsFor)
    magnification=magnificationsToExportTiffsFor{i};

    fprintf('%%%%% magnification = %s %%%%%\n',magnification)


for ii=indices
    datetime
    condName=Tcond.condName{ii}

    % get both (or one) ND2 file and export to Tiff
    %get firstRun file path
    Nd2_dir=Tcond.(['Nd2_',magnification,'_dir']){ii};

    filenametext=Tcond.(['Nd2_',magnification,'_filenametext'])(ii);

    listND2_inDir=listFilesWithExtension('.nd2','directoryPath',Nd2_dir);
    listInd=find(all([contains(listND2_inDir,filenametext),~startsWith(listND2_inDir,'.')],2));

    if length(listInd)==1
    %assert(length(listInd)==1)
    filepath=[Nd2_dir,filesep,listND2_inDir{listInd}];

    if isfile(fullfile(Nd2_dir,'channelMapInfo.xlsx'))
        useCustomChannelMap=true;
        channelMapTable=readtable(fullfile(Nd2_dir,'channelMapInfo.xlsx'));
        channelMap=containers.Map(channelMapTable{:,1}',channelMapTable{:,2}');
    else
        useCustomChannelMap=false;
    end

    % create directory for Tiffs
    outDir=Tcond.(['Tiff_',magnification,'_dir']){ii}
    if isdir(outDir)
        warning('overwriting data in %s',outDir)
        rmdir(outDir,'s')
        mkdir(outDir)
    else
        mkdir(outDir)
    end

    % output Tiffs
    if useCustomChannelMap
        nd2toTiff(filepath,'outDir',outDir,'outputXYpositions',true,'channelMap',channelMap)
    else
        nd2toTiff(filepath,'outDir',outDir,'outputXYpositions',true)
    end

    else
        warning('skipping condName=%s Mag=%s, did not find one file with filenametext in folder %s',condName,magnification,Nd2_dir)
    end
end
end
cd(parentDir)

