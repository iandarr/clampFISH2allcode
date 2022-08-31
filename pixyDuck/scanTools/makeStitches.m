function makeStitches(Scan,varargin)

% function is unfinished

% input parser
ip=inputParser;
ip.addParameter('ScanBounds','',@ischar)
ip.addParameter('SubregionArray',[1 1],@(x) and(isnumeric(x),isequal(size(x),[1 2])))
ip.addParameter('StitchSubregionSubsetOnly',[])
ip.addParameter('resolutions','native')
ip.addParameter('zplanes','all',@(x) strcmp(x,'all') || (isnumeric(x)&&isvector(x)&& all(rem(x,1)==0)) )
ip.addParameter('MaxMerge',true,@islogical)
ip.addParameter('tileType','fullTiles',@(x) ismember(x,{'fullTile','coreTile'}))
ip.addParameter('tiffOutParentDir',pwd,@isfolder)
ip.addParameter('tifNameFormat',{'R','round','_','channelNames','_','channelLabels'},@iscell)
ip.addParameter('StitchCoordsMetadataOutput',false,@islogical)
%ip.addParameter('OutputOnlyStitchCoords',false,@islogical)
ip.parse(varargin{:})

% assign variables
SubregionArray=ip.Results.SubregionArray;
StitchSubregionSubsetOnly=ip.Results.StitchSubregionSubsetOnly;
ScanBounds=ip.Results.ScanBounds;
resolutionsInput=ip.Results.resolutions;
tileType=ip.Results.tileType;
tiffOutParentDir=ip.Results.tiffOutParentDir;
tifNameFormat=ip.Results.tifNameFormat;
zplanes=ip.Results.zplanes;
MaxMerge=ip.Results.MaxMerge;
StitchCoordsMetadataOutput=ip.Results.StitchCoordsMetadataOutput;
%OutputOnlyStitchCoords=ip.Results.OutputOnlyStitchCoords;

% get ScanBounds input
ScanBoundsAllowable=fields(Scan.regions.(tileType));
if ~ischar(ScanBounds)
    error('input ScanBounds must be of type char and one of %s',join(ScanBoundsAllowable,', '))
end
if ~isempty(ScanBounds)
    if ~ismember(ScanBounds,ScanBoundsAllowable)
        error('input ScanBounds must be either: %s',join(ScanBoundsAllowable,', '))
    end
elseif ismember('innerBoxIntersect',ScanBoundsAllowable)
    ScanBounds='innerBoxIntersect';
elseif ismember('outerBoxIntersect',ScanBoundsAllowable)
    ScanBounds='outerBoxIntersect';
else
    error('specify a ScanBounds input, since default is not present in the data')
end


% get table of the stitches & their resolutions to be made
[resTable,iResStitch2iResScan]=makeResolutionTable(Scan,resolutionsInput);

% get full bounding rscan region
numResStitch=length(unique(resTable.iResStitch));
%fullReg=repmat(struct(),1,numResStitch);
for iRes=1:numResStitch
    iResScan=iResStitch2iResScan(iRes);
    fullReg(iRes).T=Scan.regions.(tileType).(ScanBounds).rscan.res(iResScan).T;
    fullReg(iRes).polyvec=Scan.regions.(tileType).(ScanBounds).rscan.res(iResScan).polyvec;
    fullReg(iRes).numRows=1+ fullReg(iRes).T.BottomRow - fullReg(iRes).T.TopRow;
    fullReg(iRes).numCols=1+ fullReg(iRes).T.RightCol - fullReg(iRes).T.LeftCol;
end

% check whether zplanes input makes sense
if ~strcmp(zplanes,'all')
for thisRound=1:Scan.numRounds
    if max(zplanes)>Scan.Rounds(thisRound).numZPerTile
        error("the zplane's max input(%i) exceeds the numZPerTile according to Scan object for round %i",max(zplanes),thisRound)
    end
end
end
%% subregions (resolution-independent information)

% array information
numSubregionArrayRows=SubregionArray(1);
numSubregionArrayCols=SubregionArray(2);
numSubregions=numSubregionArrayRows*numSubregionArrayCols;
subregionIndexMatrix=transpose(reshape(1:numSubregions,numSubregionArrayCols,numSubregionArrayRows));




% what subregions to stitch
if ~isempty(StitchSubregionSubsetOnly)
    if strcmp(StitchSubregionSubsetOnly,'corners')
        subregionIndicesToStitch=reshape(subregionIndexMatrix([1 end],[1 end])',1,4);
    elseif isnumeric(StitchSubregionSubsetOnly)
        assert(size(StitchSubregionSubsetOnly,2)==2);
        subregionIndicesToStitch=nan(1,size(StitchSubregionSubsetOnly,1));
        for i=1:size(StitchSubregionSubsetOnly,1)
            subregionIndicesToStitch(i)=subregionIndexMatrix(StitchSubregionSubsetOnly(i,1),StitchSubregionSubsetOnly(i,2));
        end
    else
        error('cannot handle this SitchSubregionSubsetOnly')
    end
else
    subregionIndicesToStitch=1:numSubregions;
end
numSubregionsToStitch=length(subregionIndicesToStitch);



fprintf('making stitches\n')
fprintf('ScanBounds for stitch will be %s\n',ScanBounds)

% numRounds
numRounds=Scan.numRounds;

% get output tif file names (resolution-independent)
sTiffFileNames=tifNameFormatInput2TifNameStruct(Scan,tifNameFormat);


% stitch
subregionCounter=0;
for subregionIndex=subregionIndicesToStitch
    subregionCounter=subregionCounter+1;
    fprintf('  starting stitch on subregionIndex=%i (%i out of %i) at %s',subregionIndex,subregionCounter,numSubregionsToStitch,char(datetime))
    if subregionCounter>1
        elapsedSeconds=toc(ticStart);
        fprintf(' (previous subregion stitch took %.0f seconds)\n',elapsedSeconds)
        ticStart=tic;
    else
        fprintf('\n');
        ticStart=tic;
    end
    
    
    [SubregionArrayRow,SubregionArrayCol]=find(subregionIndexMatrix==subregionIndex);
    
    % basic subReg structure information (all rounds)
    sAllR=struct();
    for iRes=1:numResStitch
        % size of each subregion (some bottom and right pixels of scan likely won't be used since need integer subregion sizes)
        sAllR(iRes).subregion_NumRows=floor(fullReg(iRes).numRows/numSubregionArrayRows);
        sAllR(iRes).subregion_NumCols=floor(fullReg(iRes).numCols/numSubregionArrayCols);
        
        % subregion rscan coordinates
        sAllR(iRes).rscan_subregion_TopRow =fullReg(iRes).T.TopRow + sAllR(iRes).subregion_NumRows*(SubregionArrayRow-1);
        sAllR(iRes).rscan_subregion_LeftCol=fullReg(iRes).T.LeftCol + sAllR(iRes).subregion_NumCols*(SubregionArrayCol-1);
        sAllR(iRes).rscan_subregion_BottomRow = sAllR(iRes).rscan_subregion_TopRow + sAllR(iRes).subregion_NumRows -1;
        sAllR(iRes).rscan_subregion_RightCol = sAllR(iRes).rscan_subregion_LeftCol + sAllR(iRes).subregion_NumCols -1;
        sAllR(iRes).subregionPositionVect=[sAllR(iRes).rscan_subregion_TopRow,sAllR(iRes).rscan_subregion_LeftCol,sAllR(iRes).subregion_NumRows,sAllR(iRes).subregion_NumCols];
        
        % polyshape of subregion
        sAllR(iRes).subregion_polyshape=RowColStartEnd2polyshape(sAllR(iRes).rscan_subregion_TopRow,sAllR(iRes).rscan_subregion_BottomRow,sAllR(iRes).rscan_subregion_LeftCol,sAllR(iRes).rscan_subregion_RightCol);
        assert(sAllR(iRes).subregion_polyshape.NumRegions==1);
        
    end
    
    %% channel stitching and writing to .tif

    for thisRound=1:numRounds
        s=struct(); % this subregion, and round-specific
        resTable_ThisRound=resTable(resTable.thisRound==thisRound,:);
        numResThisRound=height(resTable_ThisRound);
        
        if strcmp(zplanes,'all')
            iZ=1:Scan.Rounds(thisRound).numZPerTile;
        else
            iZ=zplanes;
        end
        if MaxMerge
            numZstitch=1;
        else
            numZstitch=length(iZ);
        end
        
        for iS=1:numResThisRound %round and res-specific coordinates loop
            iResStitch=resTable_ThisRound.iResStitch(iS);
            iResScan=resTable_ThisRound.iResScan(iS);
            
            s(iS).iResStitch=iResStitch;
            s(iS).iResScan=iResScan;
            s(iS).isNativeRes=resTable_ThisRound.isNativeRes(iS);
            
            % what are tiles overlapping this subregion
            polyshapeSubregion=sAllR(iResStitch).subregion_polyshape;
            polyvecTilesAll=Scan.Rounds(thisRound).(tileType).rscan.res(iResScan).polyvec;
            polyvecTilesOverlap=intersect(polyvecTilesAll,polyshapeSubregion);
            tileIDs=find([polyvecTilesOverlap(:).NumRegions]>0)';
            polyvecTilesOverlap=polyvecTilesOverlap(tileIDs);
            assert(all([polyvecTilesOverlap(:).NumRegions]==1))
            
            s(iS).numTilesOverlappingThisSubregion=length(tileIDs);
            
            
            % get overlapping rscan coord and subregion coord for each overlapping tile
            
            % turn polyvecTilesOverlap into table with TopRow,BottomRow,LeftCol, RightCol
            Trscan_overlapTiles=polyshape2Table(tileIDs,polyvecTilesOverlap);
            
            % make Timage_overlapTiles
            Trscan_orig=Scan.Rounds(thisRound).(tileType).rscan.res(iResScan).T;
            Timage_orig=Scan.Rounds(thisRound).(tileType).image.T;
            nativePxPerRscanPx=Scan.Rounds(thisRound).fullTiles.rscan.res(iResScan).nativePxPerRscanPx;
            Timage_overlapTiles=cropImageCoords(Trscan_orig,Trscan_overlapTiles,Timage_orig,nativePxPerRscanPx);
            
            % convert to rscan to subregion coordinates (make Tsubregion_overlapTiles)
            Tsubregion_overlapTiles=table(); Tsubregion_overlapTiles.tileID=Trscan_overlapTiles.tileID;
            Tsubregion_overlapTiles{:,{'TopRow','BottomRow'}}=round(Trscan_overlapTiles{:,{'TopRow','BottomRow'}} - sAllR(iResStitch).rscan_subregion_TopRow  + 1);
            Tsubregion_overlapTiles{:,{'LeftCol','RightCol'}}=round(Trscan_overlapTiles{:,{'LeftCol','RightCol'}} - sAllR(iResStitch).rscan_subregion_LeftCol + 1);
            Tsubregion_overlapTiles{:,{'numRows','numCols'}}=Tsubregion_overlapTiles{:,{'BottomRow','RightCol'}} - Tsubregion_overlapTiles{:,{'TopRow','LeftCol'}} + 1;
            
            % store in s
            s(iS).Tsubregion_overlapTiles=Tsubregion_overlapTiles;
            s(iS).Timage_overlapTiles=Timage_overlapTiles;
            
            assert(all(Tsubregion_overlapTiles.tileID==Timage_overlapTiles.tileID))
            
        end % end res-specific coordinates loop
        
        % get reader for this round, fill in s(iS).subregionStitch_uint16
        reader = bfGetReader();
        reader = loci.formats.Memoizer(reader,0);
        reader.setId(Scan.Rounds(thisRound).Nd2_filepath);
        
        % channel names (from scan object, these are not double checked here
        channelNames=Scan.Rounds(thisRound).channelNames;
        numChannels=length(channelNames);

        % what resolutions this round

        %iResStitch_AllThisRound=resTable.iResStitch(resTable.thisRound==thisRound)';
        %iResScan_AllThisRound=resTable.iResScan(resTable.thisRound==thisRound)';
        %isNativeRes_AllThisRound=resTable.iResScan(resTable.thisRound==thisRound)';
        
        
        %initialize subregion stitch. 3rd dimension is z, 4th dimension is channel
        for iS=1:numResThisRound
            iResStitch=s(iS).iResStitch;
            s(iS).subregionStitch_uint16=uint16(false(sAllR(iResStitch).subregion_NumRows,sAllR(iResStitch).subregion_NumCols,numZstitch,numChannels));
            if StitchCoordsMetadataOutput
                s(iS).subregionMetaImgRow_uint16=uint16(false(sAllR(iResStitch).subregion_NumRows,sAllR(iResStitch).subregion_NumCols)); % StitchCoords Metadata: ImgRow
                s(iS).subregionMetaImgCol_uint16=uint16(false(sAllR(iResStitch).subregion_NumRows,sAllR(iResStitch).subregion_NumCols)); % StitchCoords Metadata: ImgCol
                s(iS).subregionMetaImgTile_uint16=uint16(false(sAllR(iResStitch).subregion_NumRows,sAllR(iResStitch).subregion_NumCols)); % StitchCoords Metadata: ImgTile
            end
        end
        
        cameraAngleVsRefAngle=Scan.Rounds(thisRound).cameraAngleVsRefAngle;
        rowDimMustBeFlipped=Scan.Rounds(thisRound).rowDimIsFlippedVsYDim~=Scan.Rounds(Scan.referenceRound).rowDimIsFlippedVsYDim;
        
        %go through each overlapping tile
        assert(all([s(:).numTilesOverlappingThisSubregion]==s(1).numTilesOverlappingThisSubregion))
        
        for iTile=1:s(1).numTilesOverlappingThisSubregion
            tileID=s(1).Tsubregion_overlapTiles.tileID(iTile); % just use first iS since it's the same as 
           % assert(all(tileID==[s(:).Tsubregion_overlapTiles.tileID(iTile)]))
            
            for iChannel=1:numChannels
                channelName=channelNames{iChannel};
                
                for iS=1:numResThisRound
                    isNativeRes=s(iS).isNativeRes;
                    
                    % get tile vertices (in subregioncoord) overlapOnly (in this tiles)
                    subregioncoord_RowsToFill=s(iS).Tsubregion_overlapTiles{iTile,{'TopRow'}}:s(iS).Tsubregion_overlapTiles{iTile,{'BottomRow'}};
                    subregioncoord_ColsToFill=s(iS).Tsubregion_overlapTiles{iTile,{'LeftCol'}}:s(iS).Tsubregion_overlapTiles{iTile,{'RightCol'}};
                    
                    % get tile vertices (in imgcoord) overlapOnly (in this tiles)
                    imgcoord_RowsToFill=s(iS).Timage_overlapTiles{iTile,{'TopRow'}}:s(iS).Timage_overlapTiles{iTile,{'BottomRow'}};
                    imgcoord_ColsToFill=s(iS).Timage_overlapTiles{iTile,{'LeftCol'}}:s(iS).Timage_overlapTiles{iTile,{'RightCol'}};
                    
                    
                    img_uint16=getPlaneFromNd2file(reader, tileID, channelName,'closeReaderAfterGettingPlane',false,'iZ',iZ); % added z plane functionality
                    
                    % default behavior is to get max merge
                    if length(iZ)>1 && MaxMerge
                        img_uint16=max(img_uint16,[],3);
                    end
                    
                    if cameraAngleVsRefAngle~=0
                        img_uint16=rot90(img_uint16,cameraAngleVsRefAngle/90);
                    end
                    if rowDimMustBeFlipped
                        img_uint16=flipud(img_uint16);
                    end
                    
                    img_uint16_cropped=img_uint16(imgcoord_RowsToFill,imgcoord_ColsToFill,:);
                    
                    if ~isNativeRes
                        resizedRowCol=[s(iS).Tsubregion_overlapTiles.numRows(iTile),s(iS).Tsubregion_overlapTiles.numCols(iTile)];
                        img_uint16_cropped=imresize(img_uint16_cropped,resizedRowCol);
                    end
                    
                    % fill in stitched image
                    s(iS).subregionStitch_uint16(subregioncoord_RowsToFill,subregioncoord_ColsToFill,:,iChannel)=img_uint16_cropped;
                    
                    % METADATA

                    if StitchCoordsMetadataOutput
                        img_MetaImgRow_uint16=uint16(repmat([1:Scan.Rounds(thisRound).fullTileNumRows]',1,Scan.Rounds(thisRound).fullTileNumCols));
                        img_MetaImgCol_uint16=uint16(repmat([1:Scan.Rounds(thisRound).fullTileNumCols] ,Scan.Rounds(thisRound).fullTileNumRows,1));

                        if cameraAngleVsRefAngle~=0
                            img_MetaImgRow_uint16=rot90(img_MetaImgRow_uint16,cameraAngleVsRefAngle/90);
                            img_MetaImgCol_uint16=rot90(img_MetaImgCol_uint16,cameraAngleVsRefAngle/90);
                        end
                        if rowDimMustBeFlipped
                            img_MetaImgRow_uint16=flipud(img_MetaImgRow_uint16);
                            img_MetaImgCol_uint16=flipud(img_MetaImgCol_uint16);
                        end
                        img_MetaImgRow_uint16_cropped=img_MetaImgRow_uint16(imgcoord_RowsToFill,imgcoord_ColsToFill,:);
                        img_MetaImgCol_uint16_cropped=img_MetaImgCol_uint16(imgcoord_RowsToFill,imgcoord_ColsToFill,:);

                        s(iS).subregionMetaImgRow_uint16(subregioncoord_RowsToFill,subregioncoord_ColsToFill)=img_MetaImgRow_uint16_cropped;
                        s(iS).subregionMetaImgCol_uint16(subregioncoord_RowsToFill,subregioncoord_ColsToFill)=img_MetaImgCol_uint16_cropped;
                        s(iS).subregionMetaImgTile_uint16(subregioncoord_RowsToFill,subregioncoord_ColsToFill)=repmat(tileID,size(img_MetaImgRow_uint16_cropped));
                    end
                    %% If StitchCoordsMetadataOutput




                end % end res loop
                
            end % end channel loop
            
            
        end % end tile loop
        
        reader.close
        
        %% save tiffs from s(iS).subregionStitch_uint16
        %tiffOutSubregionDir=[tiffOutParentDir,filesep,Scan.scanName];
        tiffOutParentDirWithSubregionArray=[tiffOutParentDir,filesep,'SubregionArray',num2str(numSubregionArrayRows),'x',num2str(numSubregionArrayCols)];
        subregionDirStr=['Subregion_',num2str(subregionIndex),'_r',num2str(SubregionArrayRow),'_c',num2str(SubregionArrayCol)];
        tiffOutDir=[tiffOutParentDirWithSubregionArray,filesep,subregionDirStr];
        
        if ~isfolder(tiffOutDir)
            mkdir(tiffOutDir)
        end
        for iS=1:numResThisRound
            %iResStitch=s(iS).iResStitch;
            iResScan=s(iS).iResScan;
            
            for iChannel=1:numChannels
                
%                 channelName=channelNames{iChannel};
%                 channelLabel=Scan.Rounds(thisRound).channelLabels{iChannel};
%                 if isempty(channelLabel)
%                     tiffname=['R',num2str(thisRound),'_',channelName,'.tif'];
%                 else
%                     tiffname=['R',num2str(thisRound),'_',channelName,'_',channelLabel,'.tif'];
%                 end
                tiffname=[sTiffFileNames(thisRound).tifNames{iChannel},'.tif'];
                if numResStitch>1
                    tiffname=['res',num2str(iResScan),'_',tiffname];
                end
                
                
                for iZout=1:numZstitch
                    if iZout==1
                        imwrite(s(iS).subregionStitch_uint16(:,:,iZout,iChannel),fullfile(tiffOutDir,filesep,tiffname))
                    else
                        imwrite(s(iS).subregionStitch_uint16(:,:,iZout,iChannel),fullfile(tiffOutDir,filesep,tiffname),'WriteMode','append')
                    end
                end

            end % end channels
            % Metadata write
            if StitchCoordsMetadataOutput
                % decide on a prefix to designate the round (+/- resolution)
                listOfChannelTiffNames=[sTiffFileNames(thisRound).tifNames(:)];
                firstUnderscores=cell2mat(cellfun(@(x) x(1),strfind(listOfChannelTiffNames,'_'),'UniformOutput',false));
                if all(firstUnderscores(1)==firstUnderscores)
                    temp=listOfChannelTiffNames{1};
                    metadataTiffPrefix=temp(1:firstUnderscores(1));
                else
                    metadataTiffPrefix='R'&num2str(thisRound)&'_';
                end
                if numResStitch>1
                    metadataTiffPrefix=['res',num2str(iResScan),'_',metadataTiffPrefix];
                end

                imwrite(s(iS).subregionMetaImgRow_uint16,fullfile(tiffOutDir,filesep,[metadataTiffPrefix,'ImgRow.tif']));
                imwrite(s(iS).subregionMetaImgCol_uint16,fullfile(tiffOutDir,filesep,[metadataTiffPrefix,'ImgCol.tif']));
                imwrite(s(iS).subregionMetaImgTile_uint16,fullfile(tiffOutDir,filesep,[metadataTiffPrefix,'ImgTile.tif']));
            end
        
    end % end round loop for stitching channels
    
end % end subregion loop



end