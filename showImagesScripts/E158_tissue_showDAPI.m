% E158_tissue_showDAPI

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
exampleImagesDir='exampleImages/E158_tissue/';

%% get DAPI images, piece together into a single stitched image.
tiffDir='/Volumes/IAND_04/E158_tissue/scan/scan1_TumorFreshFrozen/SubregionArray5x5';

downsamplingFactor=40;
imgStruct=struct();
temp=dir(tiffDir);
subregionFolderNames=[{temp.name}]';
idRowColLookup=[[1:25]',reshape(repmat(1:5,5,1),25,1),repmat(1:5,1,5)'];
for stitchNum=1:25
    subregionDirStart=sprintf('Subregion_%s_',num2str(stitchNum));
    idx=startsWith(subregionFolderNames,subregionDirStart);
    assert(sum(idx)==1)
    subregionFolderName=subregionFolderNames{idx};
    
    dapiPath=fullfile(tiffDir,filesep,subregionFolderName,filesep,'R1_DAPI.tif');
    imgInfo=imfinfo(dapiPath);
    imgStruct(stitchNum).ArrayRow=idRowColLookup(stitchNum,2);
    imgStruct(stitchNum).ArrayCol=idRowColLookup(stitchNum,3);
    assert(~isempty(regexp(subregionFolderName,['_r',num2str(idRowColLookup(stitchNum,2)),'_'],'once')))
    assert(~isempty(regexp(subregionFolderName,['_c',num2str(idRowColLookup(stitchNum,3))],'once')))
    
    imgStruct(stitchNum).Width=imgInfo.Width;
    imgStruct(stitchNum).Height=imgInfo.Height;
    img=imread(dapiPath,'tif','PixelRegion',{[1 downsamplingFactor imgInfo.Height], [1 downsamplingFactor imgInfo.Width]});
    imgStruct(stitchNum).image=img;
    imgStruct(stitchNum).downsampledWidth=size(img,2);
    imgStruct(stitchNum).downsampledHeight=size(img,1);
end
% make stitch
assert(all([imgStruct.Width]==imgStruct(1).Width))
assert(all([imgStruct.Height]==imgStruct(1).Height))
subregionWidth=imgStruct(1).downsampledWidth;
subregionHeight=imgStruct(1).downsampledHeight;
stitchSize=[5*subregionHeight 5*subregionWidth];
dapiFullScanImg=zeros(stitchSize,'uint16');
for stitchNum=1:25
    ArrayRow=imgStruct(stitchNum).ArrayRow;
    ArrayCol=imgStruct(stitchNum).ArrayCol;
    stitchRowStart= 1 + (ArrayRow-1)*subregionHeight;
    stitchColStart= 1 + (ArrayCol-1)*subregionWidth;
    stitchRowEnd=stitchRowStart + subregionHeight - 1;
    stitchColEnd=stitchColStart + subregionWidth - 1;

    dapiFullScanImg(stitchRowStart:stitchRowEnd,stitchColStart:stitchColEnd)=...
        imgStruct(stitchNum).image;

end
dapiFullScanPath=fullfile(parentDir,filesep,exampleImagesDir,filesep,['DAPI_FullScan_downsample',num2str(downsamplingFactor),'.tif']);
imwrite(dapiFullScanImg,dapiFullScanPath)
%%
[dapiFullScanRGB,scalingMinMaxList]=mergeImgsToRGB({dapiFullScanPath},{'b'},'scalingMinMax',{[0 2000]},'convertTo','uint8');
dapiFullScanRGBPath=fullfile(parentDir,filesep,exampleImagesDir,filesep,['DAPI_FullScanRGB_downsample',num2str(downsamplingFactor),'.tif']);
imwrite(dapiFullScanRGB,dapiFullScanRGBPath)

