% E158_tissue_BothMouseHumanRegions_showImgs
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
exampleImagesDir='paper/exampleImages/E158_tissue_BothMouseHumanRegions/';

exampleInputsFile='E158_tissue_BothMouseHumanRegions_ExampleImagesInputs.xlsx';
TexAll=readtable(fullfile(parentDir,filesep,exampleImagesDir,exampleInputsFile));

%tiffParentDir='/Volumes/IAND_04/E158_tissue/scan/scan1_TumorFreshFrozen/SubregionArray5x5/';

%% Panel inputs
panelID=1;
outFile=fullfile(parentDir,exampleImagesDir,['panel',num2str(panelID)]);
S=Table2makeImgPanelsInputs(TexAll(TexAll.panelID==panelID,:));
%S.scalingInput(2,:)=repmat({{0, 97}},1,size(S.scalingInput,2));

[fig,scalingMinMax]=makeImgPanels(...
    S.imgFiles,...
    S.outColors,...
    'outFile',outFile,...
    'boundingbox',S.boundingbox,...
    'scalingMinMax',S.scalingInput,...
    'layoutMatrix',S.layoutMatrix,...
    'fontSize',20,...
    'figWidth',2048,...
    'exportFigToFile',true,...
    'inPanelLabels',S.inPanelLabels);
%% Panel inputs
panelID=2;
outFile=fullfile(parentDir,exampleImagesDir,['panel',num2str(panelID),'_DAPIonly']);
S=Table2makeImgPanelsInputs(TexAll(TexAll.panelID==panelID,:));
%S.scalingInput(2,:)=repmat({{0, 97}},1,size(S.scalingInput,2));

[fig,scalingMinMax]=makeImgPanels(...
    S.imgFiles,...
    S.outColors,...
    'outFile',outFile,...
    'boundingbox',S.boundingbox,...
    'scalingMinMax',S.scalingInput,...
    'layoutMatrix',S.layoutMatrix,...
    'fontSize',20,...
    'figWidth',2048,...
    'exportFigToFile',true,...
    'inPanelLabels',S.inPanelLabels);
