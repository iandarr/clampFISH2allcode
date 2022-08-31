function  outImgStk = getPlaneFromNd2file(ND2filenameORreader, tileID, channelName,varargin)
% outImgStk = getPlaneFromNd2file(ND2filenameORreader, tileID, channelName)
% outImgStk = getPlaneFromNd2file(ND2filenameORreader, tileID, channelName,'Name',Value)
%
% Example:
%   outImgStk = getPlaneFromNd2file('myfile.nd2', 12, 'DAPI')
%
% Example (specifying iZ and iT dimensinons)
%   outImgStk = getPlaneFromNd2file('myfile.nd2', 12, 'DAPI','iZ',2,'iT',1)
% 
% Example (getting a bounding box)
% bbox= [x y width height]
%       x       is column (starting at x=1)
%       y       is row (starting at y=1)
%       width   is number of columns
%       height  is number of rows
% 
%  myBbox = [1 1 20 30] 
%  outImg = getPlaneFromNd2file('myfile.nd2', 12, 'DAPI','bbox',myBbox)
%  % would get a 30 column x 20 row image in the top-left of the image
%
% 
% Example (specifying iZ and iT dimensinons)
%   outImgStk = getPlaneFromNd2file('myfile.nd2', 12, 'DAPI','closeReaderAfterGettingPlane',false)
%               this is quicker (but may be bad for memory? unclear)
%

ip=inputParser;
ip.addParameter('iZ',1,@(x) validateattributes(x,{'numeric'},{'vector'}));
ip.addParameter('iT',1,@(x) validateattributes(x,{'numeric'},{'scalar'}));
ip.addParameter('closeReaderAfterGettingPlane',true,@islogical)
ip.addParameter('bbox',[],@(x) and(isnumeric(x),isequal(size(x),[1 4])))
ip.parse(varargin{:});
iZvect=ip.Results.iZ;
iT=ip.Results.iT;
bbox=ip.Results.bbox;
closeReaderAfterGettingPlane=ip.Results.closeReaderAfterGettingPlane;


isValidReader = @(x) isa(x, 'loci.formats.IFormatReader') && ~isempty(x.getCurrentFile());

if isValidReader(ND2filenameORreader) % this is a reader
    %firstInputType='reader';
    reader=ND2filenameORreader;
elseif isfile(ND2filenameORreader)
    % code to leverage memoizer, though not really sure if this works
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId(ND2filenameORreader);
else
    error('could not find file %s, check that path is correct from this directory. Or, provide a valid reader of class loci.formats.ChannelSeparator',ND2filenameORreader)
end

% other irrelevant stuff
% Plane check
%isValidPlane = @(p) bfTestInRange(p,'iPlane',r.getImageCount());
% Optional tile arguments check
%isValidX = @(x) bfTestInRange(x,'x',r.getSizeX());


% set reader to the right tile (multipoint)
reader.setSeries(tileID-1)

% get the channels associated with the ND2 file
Nd2fileChannels = readND2Channels2(reader,'closeReaderAfterGettingPlane',false);

% find the channel index in the ND2 file of the requested channelName
channelIdx = find(ismember(Nd2fileChannels, channelName));
if isempty(channelIdx) % in case of wrong
    channelIdx = find(ismember(lower(Nd2fileChannels), lower(channelName)));
    if (~isempty(channelIdx))&&isscalar(channelIdx)
        warning('getPlaneFromND2file: requested channel name %s is wrong case (should have been %s). Will retreive %s instead.',channelName,Nd2fileChannels{channelIdx},Nd2fileChannels{channelIdx})
    end
end

if isempty(channelIdx)
    error('could not find channel %s in Nd2 file. Try one of the channels %s',channelName,Nd2fileChannels{:})
elseif ~(numel(channelIdx)==1)
    error('more than one channel in the Nd2 file the name %s was found',channelName)
end

%iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
% reader.getPixelType
% reader.getRGBChannelCount
for i=1:length(iZvect)
    iZ=iZvect(i);
    
    iPlane = reader.getIndex(iZ-1, channelIdx-1, iT-1) + 1;
    
    if isempty(bbox)
        planeImg  = bfGetPlane(reader, iPlane);
    else
        planeImg  = bfGetPlane(reader, iPlane, bbox(1),bbox(2),bbox(3),bbox(4));
    end
    
    if i==1
        if length(iZvect)==1
            outImgStk=planeImg;
        else
            outImgStk=zeros(size(planeImg,1),size(planeImg,2),length(iZvect),class(planeImg));
            outImgStk(:,:,i)=planeImg;
        end
    else
        outImgStk(:,:,i)=planeImg;
    end
    
end
% close reader
if closeReaderAfterGettingPlane % closing the reader takes time but it may be important for memory purposes?
    reader.close()
end

end