function channels = readND2Channels2(ND2filenameORreader,varargin)
% from dentist2 repository https://github.com/arjunrajlaboratory/dentist2
% modified to accept custom channelMap and no channelMap at all

% Dictionary of commonly used channels and their names in Elements (in
% Raj and Shaffer labs)
%    channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'GFP', 'CY3', 'Cy3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR'},...
%                            {'trans'      , 'dapi', 'gfp', 'gfp', 'tmr', 'tmr', 'alexa', 'cy', 'cy'  , 'nir', 'nir','nir'});
channelMap = containers.Map({'DAPI','YFP','CY5','CY5_1','CY5_3','CY5_2',  'CY5_5','CY5_4','Brightfield',  'CY5a', 'CY5b', 'CY5c', 'CY5d','CY5e'},...
    {'dapi','gfp','cy','cya',   'cyb',  'cyc',    'cyd',  'cye',  'trans',        'cya',  'cyb',  'cyc',  'cyd',  'cye'});

ip=inputParser;
ip.addParameter('mapChannelNames',false,@islogical)
ip.addParameter('channelMap',channelMap,@(x) strcmp(class(x),'containers.Map'))
ip.addParameter('closeReaderAfterGettingPlane',true,@islogical)
ip.parse(varargin{:});
channelMap=ip.Results.channelMap;
mapChannelNames=ip.Results.mapChannelNames;
closeReaderAfterGettingPlane=ip.Results.closeReaderAfterGettingPlane;


isValidReader = @(x) isa(x, 'loci.formats.IFormatReader') && ~isempty(x.getCurrentFile());
if isValidReader(ND2filenameORreader)
    %firstInputType='reader';
    reader=ND2filenameORreader;
elseif isfile(ND2filenameORreader)
    %firstInputType='filepath';
    
    % original code:
    %reader = bfGetReader(ND2filenameORreader);
    
    % new code to leverage memoizer, though not really sure if this works
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId(ND2filenameORreader);
    
else
    error('could not find file %s, check that path is correct from this directory. Or, provide a valid reader of class loci.formats.ChannelSeparator',ND2filenameORreader)
end


omeMeta = reader.getMetadataStore();

channels = cell(1, omeMeta.getPixelsSizeC(0).getValue());
for i = 1:numel(channels)
    c = omeMeta.getChannelName(0, i-1);
    if mapChannelNames
        channels{i} = channelMap(c.toCharArray');
    else % don't map channel names
        channels{i}=c.toCharArray';
    end
end


% close reader
if closeReaderAfterGettingPlane % closing the reader takes time but it may be important for memory purposes?
    reader.close()
end

end