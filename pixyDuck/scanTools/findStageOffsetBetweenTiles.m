function [imgcoord_OffsetsRowCol, stagecoord_OffsetsXY,stagecoord_CenterXY_ref,stagecoord_CenterXY_move,tileIDs_move]=findStageOffsetBetweenTiles(Nd2filepathOrReader_ref,Nd2filepathOrReader_move,UmPerPixel,tileIDs_ref,varargin)
%[imgcoord_RowOffset, imgcoord_ColOffset]=findStageOffsetBetweenTiles(Nd2filepathOrReader_ref,Nd2filepathOrReader_move,UmPerPixel,tileIDs_ref)
%   chooses tile ID's in 'move' file by looking for closest stage XY position
%
% [imgcoord_RowOffset, imgcoord_ColOffset]=findStageOffsetBetweenTiles(Nd2filepathOrReader_ref,Nd2filepathOrReader_move,UmPerPixel,tileIDs_ref,tileIDs_move)
%   can provide the optional name-value pair: 'tileIDs_move', tileIDs_move
%       if provided, will use those tileIDs from move file

ip=inputParser;
ip.addOptional('tileIDs_move',[],@(x) and(isnumeric(x),isequal(size(tileIDs_ref),size(x))))
ip.parse;
tileIDs_move=ip.Results.tileIDs_move;

% make column vector
if size(tileIDs_ref,1)==1
    tileIDs_ref=tileIDs_ref';
end

warning('function findControlPointOffsets is unfinished. Should not depend on RowIsDim, ColIsDim. And should account for angle. And better parse control points')

numCp=size(tileIDs_ref,1);

% get all XY stage positions, ref and move
[allStagecoord_CenterX_ref,allStagecoord_CenterY_ref]=nd2positions(Nd2filepathOrReader_ref);
[allStagecoord_CenterX_move,allStagecoord_CenterY_move]=nd2positions(Nd2filepathOrReader_move);


% convert to single 2-column XY. One row per tile (all tiles)
allStagecoord_CenterXY_ref=[allStagecoord_CenterX_ref,allStagecoord_CenterY_ref];    allStagecoord_CenterX_ref=[];  allStagecoord_CenterY_ref=[];
allStagecoord_CenterXY_move=[allStagecoord_CenterX_move,allStagecoord_CenterY_move]; allStagecoord_CenterX_move=[]; allStagecoord_CenterY_move=[];

% get reference control point XY positions
stagecoord_CenterXY_ref=allStagecoord_CenterXY_ref(tileIDs_ref,:);

% get tileIDs_move (if empty)
if isempty(tileIDs_move)
    tileIDs_move=nan(numCp,1);
    for iCp=1:numCp
        % need to find which tileID is non-reference
        cp_stagecoord_CenterXY_ref=stagecoord_CenterXY_ref(iCp,:);
        distanceList=pdist2(cp_stagecoord_CenterXY_ref,allStagecoord_CenterXY_move,'euclidean')';
        [nearestDistToMoveCp,tileID_move]=min(distanceList);
        tileIDs_move(iCp)=tileID_move;
    end
end

% get move control points XY positions
stagecoord_CenterXY_move=allStagecoord_CenterXY_move(tileIDs_move,:);


assert(isequal(size(tileIDs_ref),size(tileIDs_move)))
assert(isequal(size(stagecoord_CenterXY_ref),size(stagecoord_CenterXY_move)))

%             if isempty(str2num(cpInputDescription)) % text description
%                 if strcmpi(cpInputDescription,'UL')
%                     tileID_ref=scanMatrix(1,1);
%                 elseif strcmpi(cpInputDescription,'LL')
%                     tileID_ref=scanMatrix(end,1);
%                 elseif strcmpi(cpInputDescription,'UR')
%                     tileID_ref=scanMatrix(1,end);
%                 elseif strcmpi(cpInputDescription,'LR')
%                     tileID_ref=scanMatrix(end,end);
%                 else
%                     error('unknown cpInputDescription. In Tcond excel table it should say something like UL, LR')
%                 end
%             else %numerical tileID of reference scan has been specified
%                 tileID_ref=str2num(cpInputDescription);
%                 if ~ismember(tileID_ref,1:Scan.Rounds(roundRef).numTiles)
%                     error('cpInputDescription in Tcond has a tileID number that is not in the range of the reference scan tileIDs')
%                 end
%
%             end


imgcoord_OffsetsRowCol=nan(numCp,2);
stagecoord_OffsetsXY=nan(numCp,2);

for iCp=1:numCp
    
    tileID_ref=tileIDs_ref(iCp);
    tileID_move=tileIDs_move(iCp);
    
    % get image tiles for registration
    iZ=1;
    
    % get reference dapi image
    %fprintf('   round %i loop: reader for ref tileID=%i  from Nd2 file %s...',round,tileID_ref,Nd2filepathOrReader_ref);
    uint16_ref= getPlaneFromNd2file(Nd2filepathOrReader_ref, tileID_ref, 'DAPI','closeReaderAfterGettingPlane',false,'iZ',iZ);
    %fprintf('done.\n');
    % get move dapi image
    %fprintf('   round %i loop: reader for move tileID=%i from Nd2 file %s...',round,tileID_move,Nd2_filepath_move);
    uint16_move=getPlaneFromNd2file(Nd2filepathOrReader_move, tileID_move, 'DAPI','closeReaderAfterGettingPlane',false,'iZ',iZ);
    %fprintf('done.\n');
    
    %% convert uint16 arrays to matlab images
    assert(all(ismember({class(uint16_ref),class(uint16_move)},'uint16')))
    img_ref=mat2gray(uint16_ref,[0 65535]);
    img_move=mat2gray(uint16_move,[0 65535]);
    
    %% Find tform for this control point
    MetricThreshold=30;
    topStrongestSURFPoints=100;
    linearFractionOfImage=0.6;
    DisplayFeaturePoints=false;
    DisplayMatchingPointPairsUsedInTform=false;
    
    imgcoord_tform=estimateTform(img_move,img_ref,'MetricThreshold',MetricThreshold,'UseCentralFractionOfImage',linearFractionOfImage,'NumberOfStrongestSURFFeaturesToUse',topStrongestSURFPoints,'DisplayFeaturePoints',DisplayFeaturePoints','DisplayMatchingPointPairsUsedInTform',DisplayMatchingPointPairsUsedInTform);
    %imgcoord_tformT=imgcoord_tform.T;
    
    imgcoord_RowOffset=imgcoord_tform.T(3,2); % you would add this to the move row
    imgcoord_ColOffset=imgcoord_tform.T(3,1); % you would add this to the move row
    
    %% convert the linear shift components to a stagecoords offset (these numbers will later be added to the move stageXY positions)
    stagecoord_OffsetX=nan;
    stagecoord_OffsetY=nan;
    % Row offset
    
    %cameraAngle
    warning('needs to use camera angle still, could call stage2scancoords or something')
    
    RowIsDim='+y';
    ColIsDim='-x';
    
    if strcmpi(RowIsDim(end),'y') % img row relates to stage Y
        if strcmp(RowIsDim(1),'-') % negative
            stagecoord_OffsetY=-imgcoord_RowOffset;
        else % positive
            stagecoord_OffsetY=imgcoord_RowOffset;
        end
    elseif strcmpi(RowIsDim(end),'x') % img row relates to stage X
        if strcmp(RowIsDim(1),'-') % negative
            stagecoord_OffsetX=-imgcoord_RowOffset;
        else % positive
            stagecoord_OffsetX=imgcoord_RowOffset;
        end
    else error('RowIsDim should end in x or y with first character optionally + or -')
    end
    
    % Col offset
    if strcmpi(ColIsDim(end),'x') % img col relates to stage X
        if strcmp(ColIsDim(1),'-') % negative
            stagecoord_OffsetX=-imgcoord_ColOffset;
        else % positive
            stagecoord_OffsetX=imgcoord_ColOffset;
        end
    elseif strcmpi(ColIsDim(end),'y') % img col relates to stage Y
        if strcmp(ColIsDim(1),'-') % negative
            stagecoord_OffsetY=-imgcoord_ColOffset;
        else % positive
            stagecoord_OffsetY=imgcoord_ColOffset;
        end
    else error('ColIsDim should end in x or y with first character optionally + or -')
    end
    
    % now convert units to microns
    stagecoord_OffsetX=stagecoord_OffsetX*UmPerPixel;
    stagecoord_OffsetY=stagecoord_OffsetY*UmPerPixel;
    
    % concatenate to double column for output
    imgcoord_OffsetsRowCol(iCp,:)=[imgcoord_RowOffset,imgcoord_ColOffset];
    stagecoord_OffsetsXY(iCp,:)=[stagecoord_OffsetX,stagecoord_OffsetY];
end % end iCp loop

end