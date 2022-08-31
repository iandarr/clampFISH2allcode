% %% call maskToSegmentedDataFile.m
% 
% %% cellpose like this
% % $ conda activate cellpose
% % $ python -m cellpose --dir /Users/iandardani/Documents/_MaritSegmentFunction --img_filter dapi002         --pretrained_model nuclei --fast_mode --no_npy --diameter 90
% % $ python -m cellpose --dir /Users/iandardani/Documents/_MaritSegmentFunction --img_filter dapi003 --do_3D --pretrained_model nuclei --fast_mode --no_npy --diameter 90
% % $
% % $ python -m cellpose --dir /Users/iandardani/Documents/_MaritSegmentFunction --img_filter dapi0015 --pretrained_model nuclei --fast_mode --no_npy --diameter 90
% % dapi002 for 2D, dapi003 for 3D
% 
% %%
% cd('/Users/iandardani/Documents/_MaritSegmentFunction');
% 
% %% Various ways to call the function:
% help maskToSegmentedDataFile
% %%
% maskImgPaths={'dapi001_cp_masks.tif';'dapi002_cp_masks.tif';'dapi003_cp_masks.tif'}; % if Matlab's current directory has the files in them. If not, then give appropriate filepaths, Eg. 'data/WellA1/dapi001_cp_masks.tif',...
% maskToSegmentedDataFile(maskImgPaths,'inDir','_MaritSegmentFunction')
% %% 
% maskToSegmentedDataFile('inDir',{'_MaritSegmentFunction','_MaritSegmentFunction/Test20X'})
% %%
% maskToSegmentedDataFile(maskImgPaths{1},'inDir',{'_MaritSegmentFunction'})
% %%
% maskToSegmentedDataFile('dapi0015_cp_masks.png','inDir',{'_MaritSegmentFunction'})
% %%
% maskToSegmentedDataFile('minSegmentationPixels',600)
% 
% 
% %% Now try running rajlabimagetools functions
% improc2.segmentGUI.SegmentGUI
% %% 
% improc2.processImageObjects()
% %%  
% improc2.launchThresholdGUI()

%%
%maskToSegmentedDataFile()
%%
%improc2.processImageObjects()

%% 
improc2.launchThresholdGUI()
%%
%improc2.segmentGUI.SegmentGUI