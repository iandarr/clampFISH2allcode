function [mSumm,Tstacks,Tplanes]=getMetadata(ND2filenameORreader,varargin) 
% gets some OME metadata associated with an ND2 file.
% needs bioformats for matlab https://www.openmicroscopy.org/bio-formats/
%
% ---------------------- OUTPUTS ---------------------- 
%   mSumm is a stucture of metadata associated with the entire file,
%   including X & Y positions of each of the 1st planes in the stacks:
%
% mSumm = 
% 
%   struct with fields:
% 
%                                    numStacks: 507
%                                 numZPerStack: 1
%                             numTimesPerStack: 1
%                                      zStepUm: 0
%                            numPlanesPerStack: 10
%                                   umPerPixel: 0.3246
%                                   numImgRows: 2044
%                                   numImgCols: 2048
%                                  numChannels: 10
%                 channelInfoReferenceStackNum: 1
%                                 channelNames: {'DAPI'  'YFP'  'CY3'  'CY3_1'  'CY3_2'  'A594'  'A594_1'  'A594_2'  'CY5'  'CY5_1'}
%                       channelExposureTimesMs: [50 100 250 500 1000 250 500 1000 100 250]
%     firstZplaneWithCompleteSetOfChannelNames: 1
%                                            X: [507×1 double]
%                                            Y: [507×1 double]
%                                 channelForXY: 'DAPI'
%                                  iPlaneForXY: 1
%                               dimensionOrder: 'XYCZT'
%
% 
%  Tstacks is a table of metadata, where each row is associated with each XYposition-channel combination
%
%  Tplanes (only output if 'quickMode' is false) is a table of metadata, where each row is associated with each XYposition-channel-Zplane combination
% 
%
% OPTIONAL NAME-VALUE ARGUMENTS:
%   quickMode
%   true (default) | false
%   
%       when true, returns Tplanes=[], since it does not cycle through every plane in
%       the ND2 file. It still cycles through all planes in the first
%       stack. When false, Tplanes is a table where each row is a
%       XYposition-channel-plane combination
%
% EXAMPLE 1:
%   get file-level metadata, where each row is a XYposition-channel
%   combination. If you have Z-stacks, this metadata is associated with the
%   first plane in the stack.
%
%   [mSumm,Tstacks]=getMetadata('myND2filename.nd2')
%
% EXAMPLE 2:
%   get Tplanes data, one row for each z plane (takes longer)
%   
%   [mSumm,Tstacks,Tplanes]=getMetadata(ND2filenameORreader,'quickMode',false)
% 
ip=inputParser;
ip.addParameter('closeReaderAfterGettingPlane',true,@islogical)
ip.addParameter('quickMode',true,@islogical)
ip.addParameter('errorOnAbnormalFlag',false,@islogical)

ip.parse(varargin{:});
closeReaderAfterGettingPlane=ip.Results.closeReaderAfterGettingPlane;
quickMode=ip.Results.quickMode;
errorOnAbnormalFlag=ip.Results.errorOnAbnormalFlag;

abnormalFlag=false;

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

numStacks= omeMeta.getImageCount; % number of series (Ie. = number of multipoints = number of Z stacks)

% initialize metadataSummmary
mSumm=struct();
mSumm.numStacks=nan;
mSumm.numZPerStack=nan;
mSumm.numTimesPerStack=nan;
mSumm.zStepUm=nan;
mSumm.numPlanesPerStack=nan; % check for uniformity
mSumm.umPerPixel=nan; % check for uniformity
mSumm.numImgRows=nan; % check for uniformity
mSumm.numImgCols=nan; % check for uniformity
mSumm.numChannels=nan; % check for uniformity


% channel stuff
channelInfoRefStack=1;
mSumm.channelInfoReferenceStackNum=channelInfoRefStack;
mSumm.channelNames={''};
mSumm.channelExposureTimesMs=[];
mSumm.firstZplaneWithCompleteSetOfChannelNames='';


% XY stuff
mSumm.X=nan(numStacks,1);
mSumm.Y=nan(numStacks,1);
mSumm.channelForXY='';
mSumm.iPlaneForXY=nan;

iPlaneCounter=0; % initialize

for iStack = 1:numStacks  % loop through each series (a series is pretty much just a Z-stack. It has Z, channel, and Time dimensions)
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
    
    % get Series-level data
    reader.setSeries(iStack - 1);
    
    numImgCols = omeMeta.getPixelsSizeX(iStack - 1).getValue(); % image width, pixels
    numImgRows = omeMeta.getPixelsSizeY(iStack - 1).getValue(); % image height, pixels
    numZslices = omeMeta.getPixelsSizeZ(iStack - 1).getValue(); % number of Z slices in this stack
    numTimes=omeMeta.getPixelsSizeT(iStack - 1).getValue(); % number of times
    
    pixelColWidthUm = omeMeta.getPixelsPhysicalSizeX(iStack-1).value(ome.units.UNITS.MICROMETER).doubleValue(); %pixel width, microns
    pixelRowHeightUm = omeMeta.getPixelsPhysicalSizeY(iStack-1).value(ome.units.UNITS.MICROMETER).doubleValue(); % pixel height, microns
    if isempty(omeMeta.getPixelsPhysicalSizeZ(iStack-1))
        zStepUm= 0;
    else
        zStepUm= omeMeta.getPixelsPhysicalSizeZ(iStack-1).value(ome.units.UNITS.MICROMETER).doubleValue(); % step size in Z (microns)
    end
    
    numPlanesInStack=omeMeta.getPlaneCount(iStack-1);
    dimensionOrder=omeMeta.getPixelsDimensionOrder(iStack-1).getValue();
    numChannels=omeMeta.getChannelCount(iStack-1);
    
    
    
    % get plane-level data
    
    
    if iStack==1
        % initialize Tstacks
        Tstacks=initializeTstacks(numStacks);
        
        % initialize Tplanes, assuming for Tplanes that numPlanesThisStack will remain the same
        if quickMode==true
            Tplanes=initializeTplanes(numPlanesInStack);
        else
            Tplanes=initializeTplanes(numPlanesInStack*numStacks);
        end
        
        foundDAPIchannel=false;
    end
    
    
    % fill in Tstacks
    Tstacks{iStack,{'stack','numImgCols','numImgRows','numZslices','pixelColWidthUm','pixelRowHeightUm','zStepUm','numPlanesInStack','numChannels','numTimes'}}=[iStack numImgCols numImgRows numZslices pixelColWidthUm pixelRowHeightUm zStepUm numPlanesInStack numChannels numTimes];
    Tstacks(iStack,{'dimensionOrder'})={dimensionOrder};
    
    iChannelOfDapiOrReferenceChannel=1; % initialize to first in case we don't find DAPI channel
    iPlaneOfDapiOrReferenceChannel=1; % initialize to first in case we don't find DAPI channel
    channelNameOfReferenceChannel=char(omeMeta.getChannelName(iStack-1,iChannelOfDapiOrReferenceChannel-1));
    
    
    % fill in Tplanes (if quickMode==false, do this only for first stack)
    if or(iStack==1,quickMode==false)
        for iPlane=1:numPlanesInStack
            %iPlane = reader.getIndex(iZ - 1, iChannel -1, iTime - 1) + 1;
            
            
            iZ=                     omeMeta.getPlaneTheZ(iStack-1,iPlane-1).getValue()+1;
            iChannel=               omeMeta.getPlaneTheC(iStack-1,iPlane-1).getValue()+1;
            
            channelName=          char(omeMeta.getChannelName(iStack-1,iChannel-1)); % (iStack-1, iChannel-1)
            
            exposureTimeMs= omeMeta.getPlaneExposureTime(iStack-1,iPlane-1).value(ome.units.UNITS.MILLISECOND).doubleValue(); %convert to double
            iTime=                  omeMeta.getPlaneTheT(iStack-1,iPlane-1).getValue()+1;
            X=omeMeta.getPlanePositionX(iStack-1,iPlane-1).value().doubleValue(); %convert to double
            Y=omeMeta.getPlanePositionY(iStack-1,iPlane-1).value().doubleValue(); %convert to double
            Z=omeMeta.getPlanePositionZ(iStack-1,iPlane-1).value().doubleValue(); %convert to double
            
            % fill in Tplanes
            iPlaneCounter=iPlaneCounter+1;
            Tplanes{iPlaneCounter,{'stack','plane','iZ','iChannel','exposureTimeMs','iTime','X','Y','Z'}}=[iStack iPlane iZ iChannel exposureTimeMs iTime X Y Z];
            Tplanes(iPlaneCounter,{'channelName'})={channelName};
            
            % get iChannelOfDapiOrReferenceChannel
            if and(foundDAPIchannel==false,contains(channelName,'DAPI','IgnoreCase',true)) % we found first channel with DAPI in name
                foundDAPIchannel=true;
                iChannelOfDapiOrReferenceChannel=iChannel;
                iPlaneOfDapiOrReferenceChannel=iPlane;
                channelNameOfReferenceChannel=channelName;
            end
        end % end iPlane loop
    end     % end iPlane loop for stagement
    
    % fill in metadataSummary.X and metadataSummary.Y
    mSumm.X(iStack)=omeMeta.getPlanePositionX(iStack-1,iPlaneOfDapiOrReferenceChannel-1).value().doubleValue();
    mSumm.Y(iStack)=omeMeta.getPlanePositionY(iStack-1,iPlaneOfDapiOrReferenceChannel-1).value().doubleValue();
    
end % end iStack loop

Tplanes=Tplanes(~isnan(Tplanes.stack),:); % in case numPlanesInStack*numStacks was not the right number of rows

% close reader
if closeReaderAfterGettingPlane % closing the reader takes time but it may be important for memory purposes?
    reader.close()
end


%% fill in rest of metadataSummary

% mSumm.numStacks
mSumm.numStacks=numStacks;

% mSumm.numPlanesPerStack
if all(Tstacks.numPlanesInStack==Tstacks.numPlanesInStack(1))
    mSumm.numPlanesPerStack=Tstacks.numPlanesInStack(1);
else
    warning('numPlanesInStack are not all the same in this file'), abnormalFlag=true;
end

% mSumm.umPerPixel
if ~all([Tstacks.pixelColWidthUm==Tstacks.pixelColWidthUm(1),Tstacks.pixelRowHeightUm==Tstacks.pixelRowHeightUm(1)],2) % self consistent
    warning('not all pixelColWidthUm and/or pixelRowHeightUm are the same in each stack'), abnormalFlag=true;
elseif ~(Tstacks.pixelColWidthUm(1)==Tstacks.pixelRowHeightUm(1)) 
    warning('pixelColWidthUm is different from pixelRowHeightUm'), abnormalFlag=true; % row and column pixel sizes are not same
else
    mSumm.umPerPixel=Tstacks.pixelColWidthUm(1);
end

%  mSumm.numImgRows
if ~all(Tstacks.numImgRows==Tstacks.numImgRows(1))
    warning('numImgRows are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.numImgRows=Tstacks.numImgRows(1);
end

%  mSumm.numImgCols
if ~all(Tstacks.numImgCols==Tstacks.numImgCols(1))
    warning('numImgCols are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.numImgCols=Tstacks.numImgCols(1);
end

% mSumm.numChannels
if ~all(Tstacks.numChannels==Tstacks.numChannels(1))
    warning('numChannels are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.numChannels=Tstacks.numChannels(1);
end

% mSumm.numZPerStack
if ~all(Tstacks.numZslices==Tstacks.numZslices(1))
    warning('numZslices are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.numZPerStack=Tstacks.numZslices(1);
end

% mSumm.numTimesPerStack
if ~all(Tstacks.numTimes==Tstacks.numTimes(1))
    warning('numZPerStack are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.numTimesPerStack=Tstacks.numTimes(1);
end

% mSumm.zStepUm
if ~all(Tstacks.zStepUm==Tstacks.zStepUm(1))
    warning('zStepUm are not all the same in all stacks'), abnormalFlag=true;
else
    mSumm.zStepUm=Tstacks.zStepUm(1);
end


% mSumm.channelNames
allChannelNamesInRefStack=unique(Tplanes.channelName(Tplanes.stack==channelInfoRefStack),'stable')'; % stable preserves order

for iZ=unique(Tplanes.iZ)'
    idxTplanesRefStackAndThisiZ=all([Tplanes.stack==channelInfoRefStack,Tplanes.iZ==iZ],2);
    if all(ismember(allChannelNamesInRefStack,Tplanes.channelName(idxTplanesRefStackAndThisiZ)))
        mSumm.channelNames=Tplanes.channelName(idxTplanesRefStackAndThisiZ)';
        mSumm.firstZplaneWithCompleteSetOfChannelNames=iZ;
        mSumm.channelExposureTimesMs=Tplanes.exposureTimeMs(idxTplanesRefStackAndThisiZ)';
        break % we have found the first iZ to have the channel names
    end
    
end

if ~(mSumm.firstZplaneWithCompleteSetOfChannelNames==1)
    abnormalFlag=true;
end

% mSumm.channelForXY
mSumm.channelForXY=channelNameOfReferenceChannel;

% mSumm.dimensionOrder
if ~all(strcmp(Tstacks.dimensionOrder,Tstacks.dimensionOrder{1}))
    warning('dimension order not the same for all stacks'), abnormalFlag=true;
else
    mSumm.dimensionOrder=Tstacks.dimensionOrder{1};
end

% Tplanes is not complete when quickMode==true (it only saw the first
% stack)
if quickMode==true
   Tplanes=[]; 
end

mSumm.iPlaneForXY=iPlaneOfDapiOrReferenceChannel;

if and(errorOnAbnormalFlag,abnormalFlag)
   error('something abnormal found about this ND2 file. Probably not all stacks have the same metadata, or not all Z have all the channels. Turn errorOnAbnormalFlag=false to ignore this and still return outputs');
end
% initialization functions for Tstacks and Tplanes tables
    function TstacksInit=initializeTstacks(numRows)
        TstacksInit=array2table(nan(numRows,10),'VariableNames',{'stack','numImgCols','numImgRows','numZslices','pixelColWidthUm','pixelRowHeightUm','zStepUm','numPlanesInStack','numChannels','numTimes'});
        TstacksInit.dimensionOrder=repmat({''},numRows,1);
        
    end

    function TplanesInit=initializeTplanes(numRows)
        TplanesInit=array2table(nan(numRows,9),'VariableNames',{'stack','plane','iZ','iChannel','exposureTimeMs','iTime','X','Y','Z'});
        TplanesInit.channelName=repmat({''},numRows,1);
        TplanesInit = movevars(TplanesInit,'channelName','After','iChannel');
    end



end