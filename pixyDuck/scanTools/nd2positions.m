function [Xpositions,Ypositions]=nd2positions(ND2filenameORreader,varargin)
% needs bioformats for matlab https://www.openmicroscopy.org/bio-formats/

ip=inputParser;
ip.addParameter('closeReaderAfterGettingPlane',true,@islogical)
ip.parse(varargin{:});
closeReaderAfterGettingPlane=ip.Results.closeReaderAfterGettingPlane;


% figure out if first input is ND2 file path or a valid reader
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

omeMeta = reader.getMetadataStore();
    
% some other things from https://docs.openmicroscopy.org/bio-formats/5.9.0/developers/matlab-dev.html
% stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
% stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
% stackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
% 
% voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
% voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
% voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeYdouble = voxelSizeY.doubleValue();                                  % The numeric value represented by this object after conversion to type double
% voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeZdouble = voxelSizeZ.doubleValue();                                  % The numeric value represented by this object after conversion to type double
%
% Zstacknumb =reader.getSeriesCount()
% omeXML = char(omeMeta.dumpXML())
    
    numMultipoints= omeMeta.getImageCount % number of Z stacks
    
    Xpositions=nan(numMultipoints,1);
    Ypositions=nan(numMultipoints,1);


        for i = 1:numMultipoints % Running through # of z stacks
            
            % example getting ome metadata:
            %   First input 'iTile-1' indicates the OME 'series'
            %       For me an OME 'series' is a Z-stack with multiple
            %       channels. It can also have multiple T (time).
            %   Second input 'iPlane-1' indicates the OME 'plane'
            %       If you have nZ z-planes, nC channels, and nT times,
            %       then you should have nZ*nC*nT planes
            %   the iPlane you want can be found like this:
            %       iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
            %       where iZ, iC, and iT are 1-indexed
            %
            % then you do something like this:
            %   planePositionZ=omeMeta.getPlanePositionZ(iTile-1,iPlane-1)
            %
            %   or this:
            %
            %   planePositionZ=omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).double()
            
            planePositionX=omeMeta.getPlanePositionX(i-1,0).value().doubleValue(); %convert to double
            planePositionY=omeMeta.getPlanePositionY(i-1,0).value().doubleValue(); %convert to double
            
            Xpositions(i)=planePositionX;% look at the first plane in stack. all (from 0 to stackSizeC*stackSizeZ-1) will have same X and Y values
            Ypositions(i)=planePositionY;
            % also save stack ID

        end
% close reader
if closeReaderAfterGettingPlane % closing the reader takes time but it may be important for memory purposes?
    reader.close()
end
end

