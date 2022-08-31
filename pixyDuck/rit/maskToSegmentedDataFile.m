function maskToSegmentedDataFile(varargin)
% generates rajlabimagtools dataXXX.mat files from
% automatically-generated segmentation masks (for example, from cellpose).
% This is a replacement for the manual improc2.segmentGUI.SegmentGUI
% 
% maskToSegmentedDataFile()
%   looks for all files in the current directory with the ending _cp_masks.tif
%   OR, if they are not present, _cp_masks.png. Will then output rajlabimagetools
%   dataXXX.mat files, where the 3-digit number XXX will match the number
%   at the start of the _cp_masks.tif or _cp_masks.png file name. 
%   maskImgs are treated as a label image, where each non-zero number
%   refers to a separate segmentation. So if the images has non-zero pixel values 1,2,and 3, then
%   3 segmentations will be created. Rajlabimagetools needs a single z-plane segmentation.
%   If a mask image is a 3D stack, a Z-projection is used to get the union of all z-slices
%   (which means nearby masks can overlap once Z-projected). Currently, for .png
%   only single-plane mask inputs are supported (use .tif if you have stacks of mask inputs)
%
%
% maskToSegmentedDataFile(maskImgPaths)
%   will use provided mask image paths maskImgPaths (those mask images
%   don't have to be in the current directory)
%
%
% maskToSegmentedDataFile(__,'Name',Value)
%     Accepts the following name-value pair arguments 
% 
%         'inDir',inDir
%         character vector | cell array of character vectors
%   
%               will find all mask images (as above) in the directory inDir instead of
%               the current directory. If you supply both maskImgPaths and inDir
%               then it will make all combinations of inDir with maskImgPaths
%
% 
%         'useZplanes',useZplanes
%         numeric vector or scalar (default to all planes)
%               
%               by default, will get a max merge of the mask in all z
%               planes (if mask image is a stack). specifying useZplanes
%               will only get the max merge of useZplanes
%
%       
%          'minSegmentationPixels',minSegmentationPixels
%           numeric non-zero whole number (default 200)
%               
%               will not output any segmentation with an area (in pixels)
%               less than minSegmentationPixels
%       
% 
%           'ignoreSmallNonconnectedAreas',TF
%           logical (true or false, default is true) | 0 or 1
%           
%               sometimes there are tiny non-connected regions (of
%               segmentations found by Matlab's bwconncomp). All but the
%               largest of these regions is discarded. Rajlabimagetools
%               doesn't like these
%
%           'overWriteDataFiles',TF
%           logical (true or false, default is true) | 0 or 1
%
%               will throw an error if you try to overwrite an existing
%               dataXXX.mat file. Turn to false if you're paranoid.
%           
% Example usage, first using cellpose to create segmentation masks
% 
%   ------- In Terminal: Getting cellpose to output masks -------
%   In Terminal, if you installed cellpose with conda and need to activate
%   the conda environment with cellpose installed, then first use the
%   command:
%
%       $ conda activate cellpose
%
%   Then, still in Terminal, call cellpose. In this case I'm doing a
%   nuclear segmentation on the dapi images. The "--img_filter dapi001" part
%   ensures that cellpose only tries to segment images whose filename
%   (without the .tif extension) is 'dapi001' rather than all the other channels. You'll of course need to
%   change the input following --dir to the director of your images, and
%   choose an appropriate diameter (in pixels) for the --diameter input
%
%       $ python -m cellpose --dir /Users/iandardani/Documents/Experiment20210810 --img_filter dapi001 --pretrained_model nuclei --fast_mode --no_npy --save_tif --diameter 80
%       $ python -m cellpose --dir /Users/iandardani/Documents/Experiment20210810 --img_filter dapi002 --pretrained_model nuclei --fast_mode --no_npy --save_tif --diameter 80
%       $ python -m cellpose --dir /Users/iandardani/Documents/Experiment20210810 --img_filter dapi003 --pretrained_model nuclei --fast_mode --no_npy --save_tif --diameter 80
%
%       consider putting these commands in a .txt file and then just running
%       the .txt file from Terminal.
%
%   which will add files named:
%       dapi001_cp_masks.tif
%       dapi002_cp_masks.tif
%       dapi003_cp_masks.tif
%
%       into the folder. Or, they will have file extension .png if you didn't choose a .tif
%
%   ------- In MATLAB: Using maskToSegmentedDataFile.m to generate data0XX.mat files used in rajlabimagetools ------
%   
%   ---- Example 1 ---- 
%
%   maskImgPaths={'dapi001_cp_masks.tif';'dapi002_cp_masks.tif';'dapi003_cp_masks.tif'}; % if Matlab's current directory has the files in them. If not, then give appropriate filepaths, Eg. 'data/WellA1/dapi001_cp_masks.tif',...
%   maskToSegmentedDataFile(maskImgPaths)
%  
%   which will add these files to the same directory of the mask images:
%    data001.mat
%    data002.mat
%    data003.mat
%   which is the same output that would have been produced after running the manual improc2.segmentGUI.SegmentGUI
%
% 
%   ----- Example 2 ----  
%   first make the current directory the path with the _cp_masks.tif or
%   _cp_masks.png files. Then run:
% 
%   maskToSegmentedDataFile()
%
% 
%   ----- Example 3 ----
%   process all files in multiple directories, and change minSegmentationPixels
% 
%   inDir={'Exp1/well1','Exp1/well2','Exp2/well1','Exp2/well2'};
%   maskToSegmentedDataFile('inDir',inDir,'minSegmentationPixels',400)
%
% 
%   ----- Example 4 ----
%   process a subset of regularly-named files in multiple directories:
%
%   maskImgPaths={'dapi001_cp_masks.png';'dapi002_cp_masks.png'};
%   inDir={'Exp1/well1';'Exp1/well2';'Exp2/well1';'Exp2/well2'};
%   maskToSegmentedDataFile(maskImgPaths,'inDir',inDir)
%   
%   which will work on the following mask files:
%     {'Exp1/well1/dapi001_cp_masks.png'}
%     {'Exp1/well1/dapi002_cp_masks.png'}
%     {'Exp1/well2/dapi001_cp_masks.png'}
%     {'Exp1/well2/dapi002_cp_masks.png'}
%     {'Exp2/well1/dapi001_cp_masks.png'}
%     {'Exp2/well1/dapi002_cp_masks.png'}
%     {'Exp2/well2/dapi001_cp_masks.png'}
%     {'Exp2/well2/dapi002_cp_masks.png'}
%
%
%  ----- In rajlabimagetools (after running maskToSegmentedDataFile) ----
%  % navigate to a directory then run:
% 
%  improc2.segmentGUI.SegmentGUI  % not necessary. But if you'd like to view the results of this function, then run this. Otherwise skip.
%                                 
%  improc2.processImageObjects()
%
%  improc2.launchThresholdGUI()

%%
p=parseMaskSegmentInputs(varargin{:});
maskImgPaths=p.Results.maskImgPaths;
inDir=p.Results.inDir;
%outDir=p.Results.outDir;
useZplanesUserInput=p.Results.useZplanes;
minSegmentationPixels=p.Results.minSegmentationPixels;
ignoreSmallNonconnectedAreas=p.Results.ignoreSmallNonconnectedAreas;
overWriteDataFiles=p.Results.overWriteDataFiles;


if isempty(maskImgPaths)
    if isempty(inDir) % look in current directory for maskImgPaths
        maskImgPaths=getMaskImgPathsFromDir(pwd);
    else % look in inDir
        maskImgPaths=getMaskImgPathsFromDir(inDir);
    end
elseif ~isempty(inDir) % inDir has something in it
    % then create maskImgPaths based on inDir and maskImgPaths,
    % assuming each of the maskImgPaths is in inDir
    maskImgPaths=cellstr(maskImgPaths); % if it isnt cell already
    inDir=cellstr(inDir);% if it isnt cell already
    if size(inDir,1)==1
        inDir=inDir'; % column vector
    end
    for i=1:length(inDir) % make ending a filesep
        inDirChar=inDir{i};
        if ~strcmp(inDirChar(end),filesep), inDir{i}=[inDirChar,filesep]; end % add filesep if needed
    end
    tempMaskImgPaths=cat(3,repmat(inDir',length(maskImgPaths),1),...
        repmat(maskImgPaths,1,length(inDir)));
    tempMaskImgPaths=join(tempMaskImgPaths,'',3);
    maskImgPaths=reshape(tempMaskImgPaths,numel(tempMaskImgPaths),1);
else % just find mask images in inDir
    maskImgPaths=getMaskImgPathsFromDir(inDir);
end

% check that we can find maskImgPaths

for i=1:length(maskImgPaths)
    if ~isfile(maskImgPaths{i})
        
    end
end

numMaskFiles=length(maskImgPaths);

% if ~isempty(outDir) && ~isfolder(outDir)
%    fprintf('could not find outDir=%s, so it will be created',outDir)
%    mkdir(outDir)
% end

numDataFilesSaved=0; % initialize a count

if numMaskFiles>0
    fprintf('maskToSegmentedDataFile: finding creating rajlabimagetools masks for %i files\n',numMaskFiles)
end
for iFile=1:numMaskFiles
    fprintf('maskToSegmentedDataFile: working on file %i of %i\n',iFile,numMaskFiles)
    maskImgPath=maskImgPaths{iFile};
    if ~isfile(maskImgPath)
        warning('maskToSegmentedDataFile: could not find file %s',maskImgPath)
        continue
    end
    
    imgInfo=imfinfo(maskImgPath);
    imgFormat=imgInfo.Format; % 'tif' or 'png'
    numZ=length(imgInfo);
    if ~isempty(useZplanesUserInput)
        if ~all(ismember(useZplanesUserInput,1:numZ))
            error('user specified useZplanes=[%s] but there are only %i planes in the image %s',num2str(useZplanes),numZ,maskImgPath)
        end
        useZplanes=useZplanesUserInput;
    else
        useZplanes=1:numZ;
    end
    
    if strcmp(imgFormat,'png') && ~(isscalar(useZplanes) && useZplanes==1)
        error('can only read single-plane png files. If you need a stack of masks, then make it a tif');
    end
    
    assert(all(imgInfo(1).Width == [imgInfo.Width]));
    assert(all(imgInfo(1).Height == [imgInfo.Height]));
    imWidth=imgInfo(1).Width;
    imHeight=imgInfo(1).Height;
    
    % get image (in a stack if there is one)
    for iZ=1:length(useZplanes)
        zplane=useZplanes(iZ);
        switch imgFormat
            case 'tif'
                img1p=imread(maskImgPath,zplane);
            case 'png'
                img1p=imread(maskImgPath,'BackgroundColor','none');
            otherwise
                error('cannot handle formats that arent tif or png')
        end
        
        % make it a stack
        if iZ==1
            img=img1p;
        else
            img=cat(3,img,img1p); % not sure if img stacking them right
        end
    end
    
    % figure out how many segmentations are here
    maxSegID=max(img,[],'all');
    
    % get data00X.mat file name for output
    fileNumStr=getFileNumStr(maskImgPath);
    [maskFileDir,~,~]=fileparts(maskImgPath);
    dataOutfilePath=fullfile(maskFileDir,filesep,['data',fileNumStr,'.mat']);
    objectDirAbsPath=getAbsolutePath(maskImgPath);
    % delete existing data file if it's there
    if isfile(dataOutfilePath)
        if  overWriteDataFiles
            delete(dataOutfilePath)
        else
            error('%s already exists and overWriteDataFiles==false. Delete the file or set overWriteDataFiles=true',dataOutfilePath)
        end
    end
    
    objects=[]; % initialize objects
    % get each segmentation object
    for segID=1:maxSegID
        reportInterval=10;
        if rem(segID,reportInterval)==0
            fprintf('     working on segmentation %i of %i (reporting status every %i segmentations)\n',segID,maxSegID,reportInterval)
        end
        segMask= segID == img;
        if size(segMask,3)>1
            segMask=max(segMask,[],3); % if it is a stack, then get max merge
        end
        
        % check the mask for multiple regions
        if ignoreSmallNonconnectedAreas % don't think rajlabimagetools can handle things when there's actually 2+ regions
            segMask=removeNonconnectedAreas(segMask);
        end
        
        % check if segmentation is big enough
        segPixels=sum(segMask,'all');
        segIsBigEnough=segPixels>=minSegmentationPixels;
        if ~segIsBigEnough
            warning('discarded segmentation %i because it was only %i pixels in area (< minSegmentationPixels of %i), from file %s',segID,segPixels,minSegmentationPixels,maskImgPath)
        
        else
            % add newObj to objects
            newObj=improc2.buildImageObject(segMask,fileNumStr,objectDirAbsPath); % 
            objects=[objects,newObj];
        end
        
    end % end segmentation loop
    
    % save objects to dataXXX.mat file
    if ~isempty(objects)
    save(dataOutfilePath,'objects');
    numDataFilesSaved=numDataFilesSaved+1;
    else
        warning('did not make any rajlabimagetools objects from %s ... something different about that file',maskImgPath)
    end
end % end mask file loop



fprintf('maskToSegmentedData finished saving %i rajlabimagetools data files\n',numDataFilesSaved)
end


function p=parseMaskSegmentInputs(varargin)

p=inputParser;
p.addParameter('maskImgPaths','',@(x) ischar(x) || (iscell(x) && size(x,2)==1));
p.addParameter('inDir','',@(x) (ischar(x) || iscell(x)) && isvector(x));
%p.addParameter('outDir','',@ischar)
p.addParameter('ignoreSmallNonconnectedAreas',true,@islogical); % if the segmentation mask has a non-connecting region (an 'island'), then ignore this. Only the largest connected region for a given mask label will be used
p.addParameter('minSegmentationPixels',200,@(x) isnumeric(x) && isscalar(x))
p.addParameter('overWriteDataFiles',true,@islogical);
p.addParameter('useZplanes',[],@(x) isnumeric(x) && isvector(x));


if nargin>0
    % check over first input which may be maskImgPaths but may not follow the name-value format
    if (ischar(varargin{1}) || iscell(varargin{1})) && ~any(ismember(varargin{1},p.Parameters))
        % then maskImgPaths provided for input 1
        tempMaskImgPaths=cellstr(varargin{1});
        %         for i=1:length(tempMaskImgPaths) % check over maskImgPaths
        %             if ~isfile(tempMaskImgPaths{i})
        %                 error("first input is interpreted to be maskImgPaths. Could not find file %i %s though. Check that you're in the right directory or use 'inDir',inDir name-value pair argument",i,tempMaskImgPaths{i})
        %             end
        %         end
        % updated varargin so inputParser can understand name-value pair
        varargin=[{'maskImgPaths'},varargin];
    elseif ~any(ismember(varargin{1},p.Parameters)) % check that first input is even okay
        validInputsStr=join(p.Parameters,'  ');
        error('first input to maskToSegmentedDataFile should either be a valid maskImgPaths or a name-value from this list: %s',validInputsStr)
    end
end
p.parse(varargin{:});

end

function maskImgPaths=getMaskImgPathsFromDir(maskDirAll)

if ischar(maskDirAll)
    maskDirAll=cellstr(maskDirAll);
end
assert(isvector(maskDirAll))

numDir=length(maskDirAll);
maskImgPaths={};
for  iDir=1:numDir
    maskDir=maskDirAll{iDir};
    if ~ischar(maskDir)
        error('provided maskDir must be a character vector or a cell array of character vectors');
    elseif  ~isfolder(maskDir)
        error('could not find directory %s from current directory %s',maskDir,pwd);
    end
    temp=dir(maskDir);
    dirContents={temp.name}';
    dirContents=dirContents(~startsWith(dirContents,'.')); % no hidden files
    
    % try finding files ending in _cp_masks.tif
    idxMaskFilesTif=endsWith(dirContents,'_cp_masks.tif');
    idxMaskFilesPng=endsWith(dirContents,'_cp_masks.png');
    
    if sum(idxMaskFilesTif)>0
        maskImgPathsThisDir=dirContents(idxMaskFilesTif);
        if sum(idxMaskFilesPng)>0
            warning('maskToSegmentedDataFile is ignoreing ignoring %i files ending in _cp_masks.png in directory %s, since %i files ending in _cp_masks.tif were found',sum(idxMaskFilesPng),maskDir,sum(idxMaskFilesTif))
        end
    elseif sum(idxMaskFilesPng)>0
        maskImgPathsThisDir=dirContents(idxMaskFilesPng);
    else
        warning('maskToSegmentedDataFile did not find any images in folder %s',maskDir);
    end
    % add on directory
    maskImgPathsThisDir=fullfile(maskDir,filesep,maskImgPathsThisDir);
    
    maskImgPaths=[maskImgPaths;maskImgPathsThisDir];
end
end

function fileNumStr=getFileNumStr(maskImgPath)
[~,maskFilename,~]=fileparts(maskImgPath);
%maskFilename='dapi00065_01_cp_masks'
digitPos=regexp(maskFilename,'[0-9]');
diffQC=diff(digitPos);
if ~all(diffQC(1:2)==1)
    error('could not determine the 3-digit file number from %s. Looking for initial 3 consecutive digits',maskImgPath)
else %length(digitPos)>3 % then have 4+ digits consecutively
    if ~all(diffQC==1)
        numInitDigits=find(diffQC~=1,1);
    else
        numInitDigits=length(digitPos);
    end
end
fileNumStr=maskFilename(digitPos(1:numInitDigits));
end

function segMask=removeNonconnectedAreas(segMask)
        CC=bwconncomp(segMask);
        if CC.NumObjects>1
            fprintf('   removing nonconnected areas for a mask\n') 
            % then we have 2 or more objects
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggestNumPixels,idxBiggest] = max(numPixels);
            segMask=false(size(segMask));
            segMask(CC.PixelIdxList{idxBiggest})=true; % only keep the biggest object's mask
            
            %             [smallestNumPixels,idxSmallest] = min(numPixels);
            %             maskImgSmallest=false(size(maskImg));
            %             maskImgSmallest(CC.PixelIdxList{idxSmallest})=true;
            
        end
end


function absDir=getAbsolutePath(pathIn)
[dirIn,~,~]=fileparts(pathIn);
% if absolute it has the filepsep / in front
if strcmp(dirIn(1),filesep)
   % then its an absolute path 
   absDir=dirIn;
else % dirIn is relative path. Make it absolute.
    absDir=fullfile(pwd,filesep,dirIn);    
end
% 
% % make pathIn absolute, if it isn't already
% 
% currentFolderParts=regexp(pwd,filesep,'split');
% pathInParts=regexp(pathIn,filesep,'split');
% if any(ismember(currentFolderParts,pathInParts))
%     % then it's probably an absolute path
%     queryPathAbsolute=queryPath;
% else
%     queryPathAbsolute=fullfile(pwd,queryPath);
%     
%     if ~isfile(queryPathAbsolute)
%         error('cannot figure out how to make queryPath an absolute path, so just input it as an absolute path starting /Users..')
%     end
% end
end