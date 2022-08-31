%% call_makeScanObject


% BeforeStrip60X, AfterStrip60X,
Nd2FilepathList={...
    '/Volumes/IAND_04/E157_validation/ReadoutStrip/60X/s7_60X_BeforeStrip/20210311_234534_240__Well2_ChannelDAPI,YFP,YFP_1,YFP_2,CY3,CY3_1,CY3_2,A594,A594_1,A594_2,CY5,CY5_1,CY5_2,Brightfield_Seq0002.nd2';
    '/Volumes/IAND_04/E157_validation/ReadoutStrip/60X/s10_60X_AfterStrip/c1w2 20210323_185940_470__ChannelDAPI,YFP,YFP_1,YFP_2,CY3,CY3_1,CY3_2,A594,A594_1,A594_2,CY5,CY5_1,CY5_2,Brightfield_Seq0000.nd2';...
    '/Volumes/IAND_04/E157_validation/ReadoutStrip/20X/s6_BeforeStrip_5x5/c1 20210311_183328_436__Well2_ChannelDAPI,YFP,YFP_1,YFP_2,YFP_3,CY3,CY3_1,CY3_2,CY3_3,A594,A594_1,A594_2,A594_3,CY5,CY5_1,CY5_3,CY5_2,Brightfield_Seq0001.nd2';...
    %'/Volumes/IAND_04/E157_validation/ReadoutStrip/20X/s9_AfterStrip_5x5/c1 20210323_164312_486__Well2_ChannelDAPI,YFP,YFP_1,YFP_2,YFP_3,CY3,CY3_1,CY3_2,CY3_3,A594,A594_1,A594_2,A594_3,CY5,CY5_1,CY5_3,CY5_2,Brightfield_Seq0001.nd2';...
    '/Volumes/IAND_04/E157_validation/ReadoutStrip/20X/s8_BeforeStrip_1x5/c1 20210323_142701_441__Well2_ChannelDAPI,YFP,YFP_1,YFP_2,YFP_3,CY3,CY3_1,CY3_2,CY3_3,A594,A594_1,A594_2,A594_3,CY5,CY5_1,CY5_3,CY5_2,Brightfield_Seq0001.nd2'};

scanDir='/Volumes/IAND_04/E157_validation/ReadoutStrip/scan';

roundNames={'BeforeStrip60X';...
    'AfterStrip60X';...
    'BeforeStrip20X';...
    'AfterStrip20X'};

referenceRound=3;
registrationRoundPairs=...
    [1 2;...
     2 4];

 % featureRecognition for cameraAngle and 
 %tileIDpairsForAngleUserInput=[1 2; 5 6]; % gives me tileIDs 3 and 4 as finding angle
 %tileIDpairsForAngleUserInput=[3 4; 4 5; 12 13; 14 15; 19 20]; % from previous successfuls only, but then 3 and 4 didn't work
%tileIDpairsForAngleUserInput=[1,10;2,9;3,8;4,7;6,15;7,14;8,13;9,12;11,20;12,19;13,18;14,17;16,25;17,24;18,23;19,22]; % ACTUALLY Ydim pairs (would be rowFlip inputs)
%tileIDpairsForAngleUserInput=[4 5; 18 19; 24 25;9 10;14 15; 19 20; 21 22]; % should be good angle ones (except 9 10 is mistaken)

cameraOrientationSettings=struct();
%cameraOrientationSettings.DisplayMatchingPointPairsUsedInTform=true;
%cameraOrientationSettings.MetricThreshold=50; % default is 100 now
% cameraOrientationSettings.numTilePairsToTryToFindAngleFor=20;
% cameraOrientationSettings.ROI_img_fixed=[1 1 1024 150]; % a smaller ROI has faster feature detection
% cameraOrientationSettings.ROI_img_move=[1 1022-150 1024 150]; % a smaller ROI has faster feature detection
% cameraOrientationSettings.NumOctaves=4;
% cameraOrientationSettings.tileIDpairsForAngleUserInput=tileIDpairsForAngleUserInput;
% cameraOrientationSettings.numSuccessfulTilePairsAfterWhichToStopFindingAngle=5; % 5


% registration of control points settings
registrationSettings=struct();
%registrationSettings.DisplayMatchingPointPairsUsedInTform=true;
%registrationSettings.MetricThreshold=100; % 100 is now default here
%registrationSettings.intensityNormType='EachImage'; %'none','EachImage',or 'SameScaling'. SameScaling is now default here
%registrationSettings.featuresMatchThreshold=5; % 0<T<=100. matchFeatures internal default is 1 (for non-binary). gets passed to matchFeatures. allows T percent 'distance' in features space of a perfect match. Higher T allows more features to be matched
%registrationSettings.featuresMatchMaxRatio=0.6; %0 < R <= 1. matchFeatures internal default is 0.6. gets passed to matchFeatures. specifying a ratio threshold for rejecting ambiguous matches. Increase R to return more matches.


%%
% if exist('controlPointInitialGuessesTable','var')
%     clear controlPointInitialGuessesTable
% end
% controlPointInitialGuessesTable=table();
% controlPointInitialGuessesTable.fixedRound=[3 3 3]';
% controlPointInitialGuessesTable.moveRound =[1 4 4]';
% controlPointInitialGuessesTable.fixedX =[2 0 5]';
% controlPointInitialGuessesTable.fixedY =[2 0 5]';
% controlPointInitialGuessesTable.moveX =[0 0 10]';
% controlPointInitialGuessesTable.moveY =[0 0 10]';
%     
% Scan=makeScanObject(...
%     Nd2FilepathList,...
%     'referenceRound',3,...
%     'registrationRoundPairs',registrationRoundPairs,...
%     'cameraOrientationSettings',cameraOrientationSettings,...
%     'registrationSettings',registrationSettings,...
%     'cameraAngleOfReferenceRound',-179.8042,...
% 	'rowDimIsFlippedVsYDim',true);
% 
%     %'controlPointInitialGuessesTable',controlPointInitialGuessesTable);%,...
%     %'featureRecognitionSettingsStruct',featureRecognitionSettingsStruct);
% %% save Scan
% if ~isdir(scanDir)
%     mkdir(scanDir)
% end
% save(fullfile(scanDir,filesep,'ScanObject'),'Scan', '-v7.3')

%% stitch
tiffOutParentDir=scanDir;

load(fullfile(scanDir,filesep,'ScanObject.mat'))
makeStithes(Scan,'SubregionArray',[2 2],...
    'SitchSubregionSubsetOnly',[1 1; 1 2; 2 2],...
    'resolutions','all',...
    'tiffOutParentDir',tiffOutParentDir);

