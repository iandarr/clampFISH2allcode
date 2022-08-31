function [imgMergedOut,scalingMinMaxList]=mergeImgsToRGB(imgInputList,outColors,varargin)
% rgbMergeImgOut=mergeImgsToRGB(imgInputList,outColors,'Name',value)
% takes single-channel, single-plane images and merges them in an RGB image
%
% imgInputList is a cell array of size nImgInput x 1 with values that are either images or filepaths of images
%
% outColors is a nImgInputx1 cell array, with values of either a color
%           shortcut (Eg. 'w' for white, 'r' for red) or an RGB triplet (Eg. [1 1 1] for white, [0 0 1] for blue)
%
%
% Name-value arguments suppored are as follows:
%       mergeImgsToRGB(imgInputList,outColors,'zplanes',zplanes)
%               examples of zplanes values:
%                   a specific plane of an image (zplanes=3)
%                   max of a range of planes (zplanes=1:5)
%                   max of all planes (zplanes='all')
%
%
%       mergeImgsToRGB(imgInputList,outColors,'scalingMinMax',scalingMinMax)
%               exmples of scalingMinMax values:
%                   scalingMinMax={10 99}
%                                       for each individual channel, min is 10th
%                                       percentile and max is 99th percentile
%                   scalingMinMax={{10 99};[3000 40000]} 
%                                       for first channel: min=10th percentile and max=99th percentile
%                                       for second channel: min=3000 and max=40000
%                   scalingMinMax={[nan nan];[3000 40000]} 
%                                       for first channel: min=0, max=65535
%                                       for second channel: min=3000 and max=40000
%       mergeImgsToRGB(imgInputList,outColors,'alphaList',alphaList)
%               examples of alphaList values:
%                   [1.5; 1] % to make first of 2 images 50% brighter in relation to the first
%
% 
%       mergeImgsToRGB(imgInputList,outColors,'boundingbox',boundingbox)
%              boundingbox is a 4-element vector with the form:
%               [LeftColumn TopRow   WIDTH HEIGHT]
%
%       mergeImgsToRGB(imgInputList,outColors,'convertToUint16',convertToUint16)
%              convertToUint16=true (default) converts back to uint16. 0 maps to 0 and >=1 maps to 65535
%              convertToUint16=false keeps it as double, where original images are scaled [0 1],
%                               although when combining multiple channels which add together (as this program does) pixel values can be >1
%
%   Update 13-Oct-2021: In order to speed up performance with larger images files, if boundingbox is provided then only that portion of the tiff is read in. This means that if any scalingMinMax are provided as percentiles, these will now be the percentiles of the cropped image, not the whole image.   

%% Process inputs
p = inputParser;
%p.addRequired('imgInputList',@iscell);
%p.addRequired('outColors',@iscell);
p.addParameter('zplanes','all',@(x) any([strcmp('all',x),ismatrix(x),iscell(x)]));
p.addParameter('scalingMinMax',{},@iscell)
p.addParameter('gamma',{})
p.addParameter('alpha',{})
p.addParameter('boundingbox',[])
p.addParameter('convertTo',[],@(x) isempty(x) || ismember(lower(x),{'uint16','uint8'}))
%p.addParameter('zplanes','max',@ischar)
p.parse(varargin{:});
scalingMinMaxInput=p.Results.scalingMinMax;
zplanesInput=p.Results.zplanes;
boundingboxInput=p.Results.boundingbox;
gammaListInput=p.Results.gamma;
alphaListInput=p.Results.alpha;
convertTo=p.Results.convertTo;




%% check over imgInputList and get the images into imgOrigList
assert(iscell(imgInputList));
nImgs=length(imgInputList);

%% zplanes
% make it nImgsx1 cell array if it isn't already
if ischar(zplanesInput)
    zplanes={zplanesInput};
elseif isnumeric(zplanesInput)
    if all(rem(zplanesInput,1)==0) && size(zplanesInput,1)==1
        zplanes={zplanesInput};
    else
        error('if input zplanes is numeric it should be 1xN (row vector) of whole numbers')
    end
end

if length(zplanes)==1
    zplanes=repmat(zplanes,nImgs,1);
elseif size(zplanes,1)~=nImgs || size(zplanes,2)~=1
    error('input zplanes should be a cell array of size nImgsx1 (columb vector), but it is has size %i x %i',size(zplanes,1), size(zplanes,2))
end

%% get imgListCurrent from files (or as matrices input) and take max merge
imgListCurrent=cell(nImgs,1);

for iImg=1:nImgs
    thisImgInput=imgInputList{iImg};
    img1p=[];
    numZplanesInImage=[];
    zplanesToGet=[];
    
    if ischar(thisImgInput)
        if isfile(thisImgInput)
            % then read the image. which zplanes
            numZplanesInImage=length(imfinfo(thisImgInput));

            zplanesToGet=parseZplanes(zplanes{iImg},numZplanesInImage);
            
            for zplane=zplanesToGet
                if isempty(boundingboxInput) % no bounding box
                    if zplane==zplanesToGet(1) % first plane
                        img1p=imread(thisImgInput,zplane);
                    else
                        img1p=max(cat(3,img1p,imread(thisImgInput,zplane)),[],3); % max merge
                    end
                else % use bounding box
                    checkBoundingBoxInput(boundingboxInput)
                    PixelRegion={[boundingboxInput(2), boundingboxInput(2)+boundingboxInput(4)-1],...
                        [boundingboxInput(1), boundingboxInput(1)+boundingboxInput(3)-1]};
                    if zplane==zplanesToGet(1) % first plane
                        img1p=imread(thisImgInput,zplane,'PixelRegion',PixelRegion);
                    else
                        img1p=max(cat(3,img1p,imread(thisImgInput,zplane,'PixelRegion',PixelRegion)),[],3); % max merge
                    end
                end
            end
        else
            error('could not find file %s\n',thisImgInput)
        end
        
    elseif ismatrix(thisImgInput) % the image was provided
        if ~isempty(boundingboxInput)
            checkBoundingBoxInput(boundingboxInput)
            thisImgInput=imcrop(thisImgInput,boundingboxInput);
        end
        
        
        numZplanesInImage=size(thisImgInput,3);
        zplanesToGet=parseZplanes(zplanes{iImg},numZplanesInImage);
        
        img1p=max(thisImgInput(:,:,zplanesToGet),[],3);
        %end
    end
    
    % now we have img1p, store it
    imgListCurrent{iImg}=img1p;
end

imHeight=size(imgListCurrent{1},1);
imWidth=size(imgListCurrent{1},2);

for iImg=1:nImgs
    if ~isequal(size(imgListCurrent{iImg}),[imHeight imWidth])
        error('all images must be same size')
    end
end


%imgListCurrent=img1pList;
%clear('img1pList')

%% check over optional scaling input and make scalingMinMax
assert(iscell(scalingMinMaxInput));
scalingMinMaxList=nan(nImgs,2); %what we will fill in

if isempty(scalingMinMaxInput)
    scalingMinMaxInput=cell(nImgs,1);
    for iImg=1:nImgs
        thisImgIntmax=intmax(class(imgListCurrent{iImg})); % for uint16 be 65535
        scalingMinMaxList(iImg,1:2)=[0 thisImgIntmax]; 
    end
    
else
    if ~(size(scalingMinMaxInput,1)==nImgs)
       if  (size(scalingMinMaxInput,1)==1)&&(size(scalingMinMaxInput,2)==2)
           % then only one percentile range given for all images Ie. {10
           % 99}, copy this to all rows to apply percentiles to all images
           scalingMinMaxInput=repmat({scalingMinMaxInput},nImgs,1);
       else
           error('could not interpret scalingMinMaxInput')
       end
    end
    
    for iImg=1:nImgs
        thisScalingMinMaxInput=scalingMinMaxInput{iImg};
        if iscell(thisScalingMinMaxInput)
            % this means numbers are percentiles
            
            minPercentile=thisScalingMinMaxInput{1};
            maxPercentile=thisScalingMinMaxInput{2};
            
            if numel(minPercentile)>1
                error('invalid scalingMinMax. Is it a column vector? (use colon not comma). Or a 1x2 cell array to apply to all channels')
            end
            if numel(maxPercentile)>1
                error('invalid scalingMinMax. Is it a column vector? (use colon not comma). Or a 1x2 cell array to apply to all channels')
            end
            
            assert(minPercentile>=0 && minPercentile<=maxPercentile)
            assert(maxPercentile>=minPercentile && maxPercentile<=100)
            
            thisImg=imgListCurrent{iImg};
            
            temp=prctile(thisImg,[minPercentile,maxPercentile],'all');
            scalingMinMaxList(iImg,1:2)=[temp(1),temp(2)];
            
        elseif isnumeric(thisScalingMinMaxInput)
            assert(isequal(size(thisScalingMinMaxInput),[1 2]))
            if isnan(thisScalingMinMaxInput(1))
                thisScalingMinMaxInput(1)=0;
            end
            if isnan(thisScalingMinMaxInput(2))
                thisScalingMinMaxInput(2)=intmax(class(imgListCurrent{iImg})); % for uint16 be 65535
            end
            scalingMinMaxList(iImg,1:2)=thisScalingMinMaxInput;
        end
    end
end



%% rescale images to 0 to 1
%imgRescaledList=cell(nImgs,1);

for iImg=1:nImgs
    
    if ~isa(imgListCurrent{iImg},'uint16')
        error('only handles uint16 class')
    end
    
    minValue=scalingMinMaxList(iImg,1);
    maxValue=scalingMinMaxList(iImg,2);
    assert(all([minValue>=0,minValue<=maxValue,maxValue<=65535]))
    
    imgListCurrent{iImg}=rescale(imgListCurrent{iImg},'InputMin',minValue,'InputMax',maxValue);
end

%imgListCurrent=imgRescaledList;
%clear('imgRescaledList')

%% get gammaList and transform imgListCurrent with it if gamma~=1
gammaList=nan(nImgs,1);
if isempty(gammaListInput)
    gammaList=ones(nImgs,1);
elseif (length(gammaListInput)==nImgs)&&isnumeric(gammaListInput)
    gammaList=gammaListInput;
    if ~all(gammaList>0)
        error('all gamma provided must be >0')
    end
else
    error('provided gammaList is invalid. should be numeric 1 x nImgs')
end

if any(gammaList~=1)
    for iImg=1:nImgs
        if gammaList(iImg)~=1
            imgListCurrent{iImg} = imgListCurrent{iImg}.^gammaList(iImg); % imgListCurrent is already scaled 0 to 1
        end
    end
end


%% get alphaList, which is the relative weighting of each image. Preference is not to use this.
alphaList=nan(nImgs,1);
if isempty(alphaListInput)
    alphaList=ones(nImgs,1);
elseif (length(alphaListInput)==nImgs)&&isnumeric(alphaListInput)
    alphaList=alphaListInput;
else
    error('provided alphaList is invalid')
end

%% check over user-supplied outColors and convert to RGB triplets (outColorsTriplet)
outColorsUserSupplied=outColors;
outColors=[];
assert(length(outColorsUserSupplied)==nImgs)

colorSpecMap=containers.Map('KeyType','char','ValueType','any');
colorSpecMap('w')=      [1 1 1];
colorSpecMap('white')=  [1 1 1];
colorSpecMap('k')=      [0 0 0];
colorSpecMap('black')=  [0 0 0];
colorSpecMap('r')=      [1 0 0];
colorSpecMap('red')=    [1 0 0];
colorSpecMap('g')=      [0 1 0];
colorSpecMap('green')=  [0 1 0];
colorSpecMap('b')=      [0 0 1];
colorSpecMap('blue')=   [0 0 1];
colorSpecMap('y')=      [1 1 0];
colorSpecMap('yellow')= [1 1 0];
colorSpecMap('m')=      [1 0 1];
colorSpecMap('magenta')=[1 0 1];
colorSpecMap('c')=      [0 1 1];
colorSpecMap('cyan')=   [0 1 1];


outColorsTripletList=cell(nImgs,1);

for iImg=1:nImgs
    % look at each outColor and make it a triplet if needed
    thisOutColorUserSupplied=outColorsUserSupplied{iImg};
    if ischar(thisOutColorUserSupplied)
        if ismember(thisOutColorUserSupplied,keys(lower(colorSpecMap)))
            thisOutColorTriplet=colorSpecMap(lower(thisOutColorUserSupplied));
        else
            errorStr=[sprintf('invalid outColor in element %i must be one of the following colorSpec values\n',iImg),...
                sprintf('%s ',keys(colorSpecMap))];
            error(errorStr)
        end
    elseif isnumeric(thisOutColorUserSupplied)
        if isequal(size(thisOutColorUserSupplied),[1 3])&&...
                all(thisOutColorUserSupplied<=1)&&...
                all(thisOutColorUserSupplied>=0)
            thisOutColorTriplet=thisOutColorUserSupplied;
        else
            error('if outColor element is numeric, it should be an RGB triplet with values from 0 to 1 inclusive')
        end
    else
        error('invalid outColor input')
    end
    
    outColorsTripletList{iImg}=thisOutColorTriplet;
end




%% merge rescaled images into RGB

imgMergedRGB=zeros(imHeight,imWidth,3);

for iImg=1:nImgs
    thisImg=imgListCurrent{iImg};
    outColorsTriplet=outColorsTripletList{iImg};
    
    imgRescaledRGB=alphaList(iImg) * cat(3, thisImg*outColorsTriplet(1), thisImg*outColorsTriplet(2), thisImg*outColorsTriplet(3) );
    
    imgMergedRGB=imgMergedRGB+imgRescaledRGB; % add up all the images
    
end


% if setBoundsAtZeroAndOne
%imgMergedRGB(imgMergedRGB<0)=0; % if values fall below
%imgMergedRGB(imgMergedRGB>1)=1; % if values fall above
% end


%% convert to unit16
if isempty(convertTo)
    imgMergedRGBconverted=imgMergedRGB;
elseif strcmp(lower(convertTo),'uint16')
    imgMergedRGBconverted=im2uint16(imgMergedRGB); % maps in [0 1] to [0 65535]. All vallues >=1 in imgMergedRGB get mapped to 1.
elseif strcmp(lower(convertTo),'uint8')
    imgMergedRGBconverted=im2uint8(imgMergedRGB);
end

%% crop with boundingbox if one is provided
% % 13-Oct-2021: this is now done during imread
% 
% if ~isempty(boundingboxInput)
%     
%     % bounding box format is like this:
%     % [colStart rowStart width height], or alternatively can think of it as:
%     % [jStart   iStart   width height], or alternatively can think of it as:
%     % [left     top      width height]
%     
%     iRange=boundingboxInput(2):boundingboxInput(2)+boundingboxInput(4)-1;
%     jRange=boundingboxInput(1):boundingboxInput(1)+boundingboxInput(3)-1;
%     imgMergedOut=imgMergedRGBconverted(iRange,jRange,:);
% else
%     imgMergedOut=imgMergedRGBconverted;
% end
imgMergedOut=imgMergedRGBconverted;

%% alternate way to read images not using imread but using Tiff
% warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
% t=Tiff('YFPe001.tif','r');
% warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')
% img=read(t);
% imshow(img,[0 10000])


end

function checkBoundingBoxInput(boundingboxInput)
%% check over boundingbox
if ~isempty(boundingboxInput)
    % bounding box format is like this:
    % [colStart rowStart width height], or alternatively can think of it as:
    % [jStart   iStart   width height], or alternatively can think of it as:
    % [left     top      width height]
    
    assert(isnumeric(boundingboxInput))
    assert(isequal(size(boundingboxInput),[1 4]))
    assert(all(boundingboxInput>0))
    %assert(boundingboxInput(1)+boundingboxInput(3)-1<=imWidth)
    %assert(boundingboxInput(2)+boundingboxInput(4)-1<=imHeight)
end
end

function zplanesToGet=parseZplanes(zInp_thisIm,numZplanesInImage)
            
if strcmp(zInp_thisIm,'all')
    zplanesToGet=1:numZplanesInImage;
elseif isnumeric(zInp_thisIm)
    zplanesToGet=zInp_thisIm;
elseif ischar(zInp_thisIm)
    if contains(zInp_thisIm,':')
        zInp_thisIm_split=split(zInp_thisIm,':');
        if length(zInp_thisIm_split)~=2 || isempty(zInp_thisIm_split{1})
            error("trouble parsing z input %s. Try something like 'all' or '3:end' or 2 or '3:5'",zInp_thisIm)
        end
        
        if isempty(zInp_thisIm_split{2})
            zInp_thisIm_split{2}='end';
        end
        if strcmp(zInp_thisIm_split{2},'end')
            zInp_thisIm_split{2}=num2str(numZplanesInImage);
        end
        
        % convert to numeric
            if ~isempty(str2num(zInp_thisIm_split{1})) && ~isempty(str2num(zInp_thisIm_split{2}))
                zplanesToGet=str2num(zInp_thisIm_split{1}) : str2num(zInp_thisIm_split{2});
            else
                error("trouble parsing z input %s. Try something like 'all' or '3:end' or 2 or '3:5'",zInp_thisIm)
            end
    elseif ~isempty(str2num(zInp_thisIm))
        zplanesToGet=str2num(zInp_thisIm);
    else
        error("trouble parsing z input %s. Try something like 'all' or '3:end' or 2 or '3:5'",zInp_thisIm)
    end
    
end

% check over validity of the numeric zplanesToGet
if ~all(rem(zplanesToGet,1)==0)
    error('not all zplanes in input were whole numbers')
end

idx_validZplanes = zplanesToGet>=1 & zplanesToGet<=numZplanesInImage;
if sum(idx_validZplanes)~=length(zplanesToGet)
    warning('only %i out of the %i z planes requested are valid. There are %i z planes in the image stack',sum(idx_validZplanes),length(zplanesToGet),numZplanesInImage)
end
if sum(idx_validZplanes)==0
    error('no valid planes were requested in zplanes input = %s in image with %i planes',zInp_thisIm,numZplanesInImage)
end
zplanesToGet=zplanesToGet(idx_validZplanes);

assert(size(zplanesToGet,1)==1) % must be row vector

% if max(zplanesToGet)>numZplanesInImage
%     error("the largest zplane that we're attempting to access (%i) is greater than number of planes in the image (%i)",max(zplanesToGet),numZplanesInImage)
% end


end