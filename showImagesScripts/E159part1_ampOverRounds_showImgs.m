% E159part1_ampOverRoudns_showImgs
%
%
% Uses 
%
% get cropped + DAPI-overlaid images of 60X, 20X, and 10X at various rounds
%

warning('change pixel width to be based on Timaging table')
parentDir='/Volumes/IAND_04/';
parentDirForExport='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
cd(parentDir)

Tcond=readtable('experiments/E159part1_amp.xlsx')
Timaging=readtable('experiments/E159part1_amp_ImagingDetails.xlsx')
exampleImagesDir='exampleImages/amplification/';
TroiIn=readtable(fullfile(exampleImagesDir,filesep,'amplification_exampleImagesInputs.xlsx'));
cd(parentDir)

drawHigherMagFOVsInsideLowerMagWholeFOV=true;
%condList=Tcond.condID;
condList=1;

outImgDir='exampleImages/amplification/';

% TroiIn --> TroiFinal
% TroiFinal is just like TroiIn but additionally has data in these columns:
%   pixelSizeMicrons
%   tiffDir
%
%   ch1ExposureTimeMs
%   ch2ExposureTimeMs
%   ... (for all channels) ...
%   ch1TiffFilepath
%   ch2TiffFilepath
%   ... (for all channels) ...
%
% which are found by referencing Tcond and Timaging
%
% In addition, if ROIs are not fully defined with image & boundingbox, then look at these columns:
%   roiNumReference
%   offsetFromReference_x
%   offsetFromReference_y
% and add that to TroiFinal

assert(all(Tcond.pixelSize_10X==Tcond.pixelSize_10X(1)))
assert(all(Tcond.pixelSize_20X==Tcond.pixelSize_20X(1)))
assert(all(Tcond.pixelSize_60X==Tcond.pixelSize_60X(1)))


assert(isequal(TroiIn.roiNum,[1:height(TroiIn)]'));
if any([TroiIn.roiNumReference>[1:height(TroiIn)]'],1)
    error('If a roi references another one, it should be lower in the table (TroiIn) than the reference roi')
end
TroiFinal=TroiIn;

for i=1:height(TroiIn)
    
    assert(strcmp(TroiIn.experiment{i},'E159'))
    condID=TroiIn.condID(i);
    rowInTcond=find(Tcond.condID==condID);
    assert(isscalar(rowInTcond));
    
    
    numChannels=TroiIn.numChannels(i);
    
    magnification=TroiIn.magnification{i};
    
    tiffDir=Tcond.(['Tiff_',magnification,'_dir']){rowInTcond};
    TroiFinal.tiffDir{i}=tiffDir;
    
    pixelSizeMicrons=Tcond.(['pixelSize_',magnification])(rowInTcond);
    
    % may need to map a reference roi to this roi
    roiNumReference=TroiIn.roiNumReference(i);
    if isnan(roiNumReference)
        % then ROI is reference --> we can get these roi filepaths now
        imgFileNum=TroiIn.imgFileNum(i);
        boundingBox=[TroiIn.bbox1(i),TroiIn.bbox2(i),TroiIn.bbox3(i),TroiIn.bbox4(i)];
        
    elseif and(isnumeric(roiNumReference),ismember(roiNumReference,TroiIn.roiNum))
        % then not reference --> this ROI must be defined by another ROI
        
        % first get information on reference ('source') roi
        % s= source or reference of ROI (Ie. a 60X image number and specified bounding box)
        % d= destination magnifications (Ie. 20X or 10X) in which to find the same ROI
        
        iRef=roiNumReference; %roiNum = row
        
        sImgFileNum=TroiFinal.imgFileNum(iRef);
        sourceMagnification=TroiFinal.magnification{iRef};
        sPixelSizeMicrons=Tcond.(['pixelSize_',sourceMagnification])(rowInTcond);
        
        sImgSize=[TroiFinal.imgHeight(iRef),TroiFinal.imgWidth(iRef)];
        dImgSize=[TroiFinal.imgHeight(i),TroiFinal.imgWidth(i)];
        dPixelSizeMicrons=Tcond.(['pixelSize_',magnification])(rowInTcond);
        
        sBoundingBox=[TroiFinal.bbox1(iRef),TroiFinal.bbox2(iRef),TroiFinal.bbox3(iRef),TroiFinal.bbox4(iRef)];
        if any(isnan(sBoundingBox))
            error('bbox1, bbox2, box3, bbox4 should be numbers')
        end
        
        sourceDir=TroiFinal.tiffDir{iRef};
        sMasksCropped={false(sBoundingBox(4),sBoundingBox(3))}; % dummy variable
        s_tempXYdata=readtable(fullfile([sourceDir,filesep,'XYdata.csv']));
        sXimg=s_tempXYdata.Xposition(sImgFileNum);
        sYimg=s_tempXYdata.Yposition(sImgFileNum);
        
        destinationDir=TroiFinal.tiffDir{i};
        d_tempXYdata=readtable(fullfile([destinationDir,filesep,'XYdata.csv']));
        dXimgs=[d_tempXYdata.Xposition];
        dYimgs=[d_tempXYdata.Yposition];
        
        XoffsetMicron=TroiFinal.offsetFromReference_x(i);
        YoffsetMicron=TroiFinal.offsetFromReference_y(i);
        
        RowIsDim='+y';
        ColIsDim='-x';
        % call function to map source ROI onto destination magnificationimage
        [imgFileNum,boundingBox,~,~]=mapSegsToNewImgs(...
            1,sBoundingBox,sMasksCropped,sXimg,sYimg,dXimgs,dYimgs,sImgSize,dImgSize,sPixelSizeMicrons,dPixelSizeMicrons,XoffsetMicron,YoffsetMicron,RowIsDim,ColIsDim);
        
        assert(isscalar(imgFileNum))
        assert(isequal(size(boundingBox),[1 4]))
        if any(isnan(imgFileNum))
            i
            Tcond.condName{rowInTcond}
            error('at least one element of dImageFileSuffix is nan. Probably means destination roi cannot fit into any one image. Try changing offsets in TroiIn')
        end
        
        TroiFinal.imgFileNum(i)=imgFileNum;
        TroiFinal{i,{'bbox1','bbox2','bbox3','bbox4'}}=boundingBox;
        
    else
        error('roiNumReference in excel sheet TroiIn should be blank or a number')
    end
    
    
    for iChannel=1:numChannels
        % fill in TiffFilepath and ExposureTimeMs
        chStr=['ch',num2str(iChannel)];
        chXName=TroiIn.([chStr,'Name']){i};
        
        filename=[chXName,num2str(imgFileNum,'%03i'),'.tif'];
        TroiFinal.([chStr,'TiffFilepath']){i}=[Tcond.(['Tiff_',magnification,'_dir']){rowInTcond},filesep,filename]; %Eg. 'dapi'
        
        rowTimaging=find(all([strcmp(Timaging.magnification,magnification),strcmp(Timaging.RajlabimagetoolsChannelName,chXName)],2));
        assert(length(rowTimaging)==1)
        
        TroiFinal.([chStr,'ExposureTimeMs'])(i)=Timaging.ExposureTimeMs(rowTimaging);
        
    end
    
    
end

writetable(TroiFinal,fullfile(exampleImagesDir,filesep,'amplification_exampleImagesFinal.xlsx'));

%TroiFinal([5 11],{'imgFileNum','bbox1','bbox2','bbox3','bbox4'})
% merge files by calling mergeImgsToRGB, then write to exampleImagesDir
cd(parentDir)
Troi=readtable(fullfile(exampleImagesDir,filesep,'amplification_exampleImagesFinal.xlsx'));

if length(unique(Troi.roiName))~=height(Troi)
    error('roiName must be unique since this will be the output image name')
end

for i=1:height(Troi)

    roiNum=i;
    condID=TroiIn.condID(i);
    boundingbox=Troi{i,{'bbox1','bbox2','bbox3','bbox4'}};
    zplanes=Troi.zplanes(i);
    if iscell(zplanes)
        zplanes=zplanes{1};
    end
    if ~isempty(str2num(zplanes))
        zplanes=str2num(zplanes);
    end
    
    numChannels=Troi.numChannels(i);
    
    
    %channel-specific inputs
    imgInputList=cell(numChannels,1);
    outColors=cell(numChannels,1);
    alphaList=nan(numChannels,1);
    scalingMinMax=cell(numChannels,1);
    for iChannel=1:numChannels
        chStr=['ch',num2str(iChannel)];
        
        % imgInputList
        imgInputList{iChannel}=Troi.([chStr,'TiffFilepath']){i};
        
        %outColors
        outColor=Troi.([chStr,'Color']){i};
        if ~isempty(str2num(outColor))
            outColor=str2num(outColor);
        end
        outColors(iChannel)={outColor};
        
        %scalingMinMax
        scalingMinMax(iChannel)={[Troi.([chStr,'ContrastMin'])(i),Troi.([chStr,'ContrastMax'])(i)]};
        
        %alphaList
        alphaList(iChannel)=Troi.([chStr,'Alpha'])(i);
        
    end
    imgMerged=mergeImgsToRGB(imgInputList,outColors,'boundingbox',boundingbox,'zplanes',zplanes,'scalingMinMax',scalingMinMax,'alphaList',alphaList);
    
    % write merged image
    mergedImageName=[Troi.roiName{i},'.tif'];
    
    cd(parentDirForExport)
    if ~isfolder(outImgDir)
        mkdir(outImgDir)
    end
    imwrite(imgMerged,fullfile(outImgDir,filesep,mergedImageName))
    
    cd(parentDir)
end

%% 

