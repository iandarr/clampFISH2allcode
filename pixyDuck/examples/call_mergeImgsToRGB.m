
imgInputList={'gfp001.tif'; 'dapi001.tif'};
outColors={'g';'b'};
%outColors={[1 1 1],[0 0 1]};
scalingMinMax={{10 99.8};{10 99.8}}; % percentile values in {} brackets
%scalingMinMax={[1700        10000];[1600        8000]}; % actual pixel values in [ ] brackets
%scalingMinMax={{10 99.8};[1600        8000]}; % first image percentiles, second image actual values 

%zplanes=1:2; % only get max merge of first two image planes (same for both image stacks provided)
[imgMerge,contrastLims]=mergeImgsToRGB(imgInputList,outColors,'scalingMinMax',scalingMinMax);%,'zplanes',zplanes);

contrastLims
%%
imshow(imgMerge)
imwrite(imgMerge,'Merge_gfp_dapi.tif')

