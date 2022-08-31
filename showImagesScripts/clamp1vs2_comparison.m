% clamp1vs2_comparison
%%
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
clamp2_dataDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';

clamp1_dataDir='rawData/clampFISH1.0/20180206_comparisonofamp_rep1/ClampFISH_Cy5_R6_GFPposcells/20X';
clamp2_dataDir='rawData/E154_pooled_amp/Tiff/s4 20X p2/ser15_pool_1to125';

clamp1_GFP_iCy5=imread(fullfile(parentDir,filesep,clamp1_dataDir,filesep,'Cy5_3000ms_2.tif')); % 3s exposure
clamp1_DAPI=imread(fullfile(parentDir,filesep,clamp1_dataDir,filesep,'DAPI_3000ms_2.tif'));

clamp2_Z=2;
clamp2_GFP_Atto647N=imread(fullfile(parentDir,filesep,clamp2_dataDir,filesep,'cyc006.tif'),clamp2_Z); % cyc=1s exposure
clamp2_DAPI=imread(fullfile(parentDir,filesep,clamp2_dataDir,filesep,'dapi006.tif'),clamp2_Z);

exampleImagesDir='paper/exampleImages/clampFISH1_vs_clampFISH2_comparison/';

%% clampFISH 1.0 merge
clamp1_sideLen=300; % for clamp1. 2x for clamp2 since resolution doubled (13um/20 for clampFISH1.0, 6.5um/20 for clampFISH 2.0)
ImgPt=[701 201];
clamp1_cropBox=[ImgPt(1) ImgPt(2) clamp1_sideLen  clamp1_sideLen];% [XMIN YMIN WIDTH HEIGHT]

[clamp1_imgMerged,CLims]=mergeImgsToRGB({clamp1_GFP_iCy5;clamp1_DAPI},{'w';'b'},'scalingMinMax',{[750 2000];{10 99.9}});
clamp1_crop=imcrop(clamp1_imgMerged,clamp1_cropBox);
imwrite(clamp1_imgMerged,fullfile(parentDir,filesep,exampleImagesDir,filesep,'GFP_clampFISH1_iCy5_Round6_20X_3000msExp_650nmRes_CCD.tif'))
imwrite(clamp1_crop,     fullfile(parentDir,filesep,exampleImagesDir,filesep,sprintf('GFP_clampFISH1_iCy5_Round6_20X_3000msExp_650nmRes_CCD_%spx_cropped.tif',num2str(clamp1_sideLen))))

% clampFISH 2.0 merge
clamp2_sideLen=clamp1_sideLen*2;
ImgPt=[1000 225];
clamp2_cropBox=[ImgPt(1) ImgPt(2) clamp2_sideLen  clamp2_sideLen];

[clamp2_imgMerged,CLims]=mergeImgsToRGB({clamp2_GFP_Atto647N;clamp2_DAPI},{'w';'b'},'scalingMinMax',{[750 3000];{10 99.9}});
clamp2_crop=imcrop(clamp2_imgMerged,clamp2_cropBox);

imwrite(clamp2_imgMerged,fullfile(parentDir,filesep,exampleImagesDir,filesep,'GFP_clampFISH2_Atto647N_Round8_20X_1000msExp_325nmRes_cCMOS.tif'))
imwrite(clamp2_crop,     fullfile(parentDir,filesep,exampleImagesDir,filesep,sprintf('GFP_clampFISH2_20X_Atto647N_Round8_1000msExp_325nmRes_cCMOS_%spx_cropped.tif',num2str(clamp2_sideLen))))

% 
% figure(1)
% imshow(clamp1_imgMerged)
% figure(2)
% imshow(clamp1_crop)
% figure(3)
% imshow(clamp2_imgMerged)
% figure(4)
% imshow(clamp2_crop)
