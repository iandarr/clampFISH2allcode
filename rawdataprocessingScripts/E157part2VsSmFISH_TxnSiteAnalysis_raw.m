% E157part2VsSmFISH_TxnSiteAnalysis
% Have smFISH at 60X (Cy3) and clampFISH 2.0 at 20X (Atto647N)
% register and stitch the data. Keep native resolutions.
%
% script is based on E157part1rep2_strip_raw

%% Load Tcond
%clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
outScanParentPath=fullfile(parentDir,filesep,'rawData/E157_validation/VsSmFISH/stitches');
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E157_validation_Tcond.xlsx');
if ~all(Tcond.condID==[1:height(Tcond)]')
    error('condID must be equivalent its row number in this code')
end
Tcond=readtable(fullfile(parentDir,'paper/experiments/E157_validation_Tcond.xlsx'));
%%

%condIDlist=[12];
%condIDlist=[14]
%condIDlist=[12;14;15];
condIDlist=[16]

% make ScanObject.mat file with all rounds, including the XY positions of each tile
for iCondID=1:length(condIDlist)
    condID=condIDlist(iCondID);
    rowInTcond=find(Tcond.condID==condID);
    assert(numel(rowInTcond)==1)
    geneName=Tcond.ReadoutSetName{rowInTcond};
    condName=Tcond.condName{rowInTcond};
    scanName=[Tcond.wellID{rowInTcond},'_',condName];
    
    % make Nd2FilepathList
    refRound=1;
    
    numRounds=2;
    Nd2FilepathList=cell(numRounds,1);
    % Imaging cycle 1 (called 'round 1' by makeScanObject) is 20X clampFISH 2.0
    thisRound=1;
    mag='20X';
    [Nd2_filepath,Nd2_filename]=dirAndFilenametext2filepath(Tcond.(['Nd2_',mag,'_dir']){rowInTcond},Tcond.(['Nd2_',mag,'_filenametext']){rowInTcond});
    Nd2FilepathList{thisRound}=Nd2_filepath;
    % Imaging cycle 2 (called 'round 2' by makeScanObject) is 60X smFISH
    thisRound=2;
    mag='60X';
    [Nd2_filepath,Nd2_filename]=dirAndFilenametext2filepath(Tcond.(['Nd2_',mag,'_dir']){rowInTcond},Tcond.(['Nd2_',mag,'_filenametext']){rowInTcond});
    Nd2FilepathList{thisRound}=Nd2_filepath;    
    
    %channelLabels=repmat({{''}},numRounds,1);
    Tchan=table();
    Tchan.round=        [   1         2         2   ]';
    Tchan.channel=      {'CY5',   'CY3',   'CY5'}';
    Tchan.channel=      {'CY5',   'CY3',   'CY5'}';
    Tchan.channelLabels={[geneName,'_20Xclamp2'],[geneName,'_60XsmFISH'],[geneName,'_60Xclamp2stripped']}';
    %{'round','channel','channelLabels','channelPrefixes'}
    
    % call makeScanObject    
    Scan=makeScanObject(Nd2FilepathList,...
       'channelLabels',Tchan,...
       'referenceRound',refRound,...
       'skipGetInnerBox',true);%,...
       %'rowDimIsFlippedVsYDim',true,...
       %'cameraAngleOfReferenceRound',-179.7526);%from wellA1, refRound=1, 8pairs: -179.7526 although was using -179.792 for a while % take out 'cameraAngle' pre-specification later on
    
    Scan.scanName=scanName;
    
    scanOutDir=fullfile(outScanParentPath,filesep,scanName);
    if ~isfolder(scanOutDir)
        mkdir(scanOutDir)
    end
	save(fullfile(scanOutDir,filesep,'ScanObject.mat'),'Scan', '-v7.3')
        
end % end scans loop

% makeStitches

%StitchSubregionSu bsetList=[1 1;3 10];
%warning('subregion subset only');
%[a,b]=ind2sub([10,10],1:100);
%SubregionsAll=[b',a'];

%StitchSubregionSubsetList=SubregionsAll(1:16,:);
%StitchSubregionSubsetList=SubregionsAll(17:end,:);

%clear Scan
for iCondID=1:length(condIDlist)
    condID=condIDlist(iCondID);
    rowInTcond=find(Tcond.condID==condID);
    assert(numel(rowInTcond)==1)
    geneName=Tcond.ReadoutSetName{rowInTcond};
    condName=Tcond.condName{rowInTcond};
    scanName=[Tcond.wellID{rowInTcond},'_',condName];
    fprintf('%s\n',repmat('-',1,200))
    fprintf('%s starting makeStitches for condID=%i at %s  %s\n',repmat('-',1,80),condID,char(datetime),repmat('-',1,80))
    
    scanOutDir=fullfile(outScanParentPath,filesep,scanName);
    load(fullfile(scanOutDir,filesep,'ScanObject.mat'))
    
    %tifNameFormat={'channelPrefixes','_','channelExposureTimesMs','ms_','channelLabels'};
    tifNameFormat={'cycle','round','_','channelNames','_','channelExposureTimesMs','ms_','channelLabels'};
    SubregionArray=[2 2];
    
    makeStitches(Scan,...
        'SubregionArray',SubregionArray,...
        'tifNameFormat',tifNameFormat,...
        'tiffOutParentDir',scanOutDir,...
        'ScanBounds','outerBoxIntersect',...
        'resolutions','all',...
        'MaxMerge',false)%,...
    %         'StitchSubregionSubsetOnly',[1 2; 2 1; 2 2]
    % 
    % 'resolutions','all' (instead of default 'native') means that both
    % imaging cycles (cycle 1 20X, cycle 2 60X will also be up- or
    % down-sampled to match the other one). Only using this to view the
    % alignment.
    %'StitchSubregionSubsetOnly', StitchSubregionSubsetList);
    
end


