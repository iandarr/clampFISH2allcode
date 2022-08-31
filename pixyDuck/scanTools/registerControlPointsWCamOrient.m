function    T=registerControlPointsWCamOrient(Scan,fixedRound,moveRound,cpTilePairsList,rstageXY_fixedAllTile,stageXY_moveAllTile,cameraAngleVsRefAngle_Move,rowDimIsFlippedVsYDim_Move,registrationSettingsStruct)
%% Assumptions about camera orientation

% these control how we convert RowCol offsets to stage coordinates 
cameraAngle_Ref=Scan.cameraAngleOfReferenceRound;
rowDimIsFlippedVsYDim_Ref=Scan.Rounds(Scan.referenceRound).rowDimIsFlippedVsYDim;

% these control whether we flip & rotate the image before registering it
cameraAngleVsRefAngle_Fixed=Scan.Rounds(fixedRound).cameraAngleVsRefAngle; % either 0,90,180,270 cannot be in between
%cameraAngleVsRefAngle_Move=Scan.Rounds(moveRound).cameraAngleVsRefAngle; % either 0,90,180,270 cannot be in between

assert(rem(cameraAngleVsRefAngle_Fixed,90)==0) % must be multiple of 90
assert(rem(cameraAngleVsRefAngle_Move,90)==0)  % must be multiple of 90

rowDimIsFlippedVsYDim_Fixed=Scan.Rounds(fixedRound).rowDimIsFlippedVsYDim;
%rowDimIsFlippedVsYDim_Move=Scan.Rounds(moveRound).rowDimIsFlippedVsYDim;

%% Register control points in cpTilePairsList
transformTypeList={'rigid';'similarity';'affine'};
T=table();
iT=0; % initialize table row index

numCpPairsToFind=size(cpTilePairsList,1);

numRowsT=numCpPairsToFind*length(transformTypeList);
T.fixedRound=repmat(fixedRound,numRowsT,1);
T.moveRound=repmat(moveRound,numRowsT,1);

%get bioformats readers
reader_fixed=loci.formats.Memoizer(bfGetReader(),0);
reader_fixed.setId(Scan.Rounds(fixedRound).Nd2_filepath);
reader_move=loci.formats.Memoizer(bfGetReader(),0);
reader_move.setId(Scan.Rounds(moveRound).Nd2_filepath);

% determine if downsampling is needed
if Scan.Rounds(fixedRound).umPerPixel==Scan.Rounds(moveRound).umPerPixel
    % keep at these resolutions
    downsampleFactor_fixed=1;
    downsampleFactor_move=1;
    umPerPixelOfRegImgs=Scan.Rounds(fixedRound).umPerPixel;
elseif Scan.Rounds(fixedRound).umPerPixel<Scan.Rounds(moveRound).umPerPixel % fixed is finer resolution, so downsample it
    downsampleFactor_fixed=Scan.Rounds(fixedRound).umPerPixel/Scan.Rounds(moveRound).umPerPixel; %factor <1
    downsampleFactor_move=1;
    umPerPixelOfRegImgs=Scan.Rounds(moveRound).umPerPixel;
else % move is finer resolution, so downsample it
    downsampleFactor_fixed=1;
    downsampleFactor_move=Scan.Rounds(moveRound).umPerPixel/Scan.Rounds(fixedRound).umPerPixel; %factor <1
    umPerPixelOfRegImgs=Scan.Rounds(fixedRound).umPerPixel;
end


        
        %% cycle through each control point pair (CpPair) to find offset
        
    for iCpPair=1:numCpPairsToFind
        
        tileID_fixed=cpTilePairsList(iCpPair,1);
        tileID_move= cpTilePairsList(iCpPair,2);
        fprintf('    iCpPair %2i of %2i (fixedRound tile %2i, moveRound tile %2i)...',iCpPair,numCpPairsToFind,tileID_fixed,tileID_move)
        
        iZ_forRegistration_fixed=ceil(0.5* Scan.Rounds(fixedRound).numZPerTile); % take the middle slice of the Z stack
        iZ_forRegistration_move= ceil(0.5* Scan.Rounds(moveRound).numZPerTile); % take the middle slice of the Z stack
        
        uint16_imgFixed= getPlaneFromNd2file(reader_fixed, tileID_fixed, 'DAPI','closeReaderAfterGettingPlane',false,'iZ',iZ_forRegistration_fixed);
        uint16_imgMove=getPlaneFromNd2file(reader_move, tileID_move, 'DAPI','closeReaderAfterGettingPlane',false,'iZ',iZ_forRegistration_move);
        
        % convert uint16 arrays to matlab images
        imgFixed=mat2gray(uint16_imgFixed,[0 65535]);
        imgMove=mat2gray(uint16_imgMove,[0 65535]);
        
        
        rstageXY_fixed=rstageXY_fixedAllTile(tileID_fixed,:);
        stageXY_move=stageXY_moveAllTile(tileID_move,:);
        % SIMULATION of rotation and simulatedRowDimIsFlippedVsYDim
        % add rotation/flipping to test this algorithm
%         simulatedImagesOn=false;
%         if simulatedImagesOn
            %simulatedExtraCameraAngle=87;
            %warning('SIMULATION of image rotation or flipping is on')
            %imgMove=flip(imgMove,1); % 1=row dimension reversed, 2=col dimension reversed
            %imgFixed=imrotate(imgFixed,-simulatedExtraCameraAngle,'crop');
                
            %imgMove=imrotate(imgMove,-simulatedExtraCameraAngle,'crop');
            %         imshowpair(imgFixed_beforeTransform,img1,'Scaling','joint');
            %         title('before (green), after (purple)')
%         end
        
        
        % if images are different sizes, downsample the finer-resolution one
        if downsampleFactor_fixed~=1
            imgFixed=imresize(imgFixed,downsampleFactor_fixed);
        end
        if downsampleFactor_move~=1
            imgMove=imresize(imgMove,downsampleFactor_move);
        end
        
        % if cameraAngle or rowDimIsFlippedVsYDim different than reference
        % round, need to flip & rotate images 
        
        % flip images
        if rowDimIsFlippedVsYDim_Ref~=rowDimIsFlippedVsYDim_Fixed
            imgFixed=flip(imgFixed,1);
        end
        if rowDimIsFlippedVsYDim_Ref~=rowDimIsFlippedVsYDim_Move
            imgMove=flip(imgMove,1);
        end
        % rotate images (only in increments of 90 degrees)
        if cameraAngleVsRefAngle_Fixed~=0
            imgFixed=rot90(imgFixed,cameraAngleVsRefAngle_Fixed/90); % 1 rotates 90degrees counterclockwise,
        end
        if cameraAngleVsRefAngle_Move~=0
            imgMove=rot90(imgMove,cameraAngleVsRefAngle_Move/90); % 1 rotates 90degrees counterclockwise,
        end
        %% call estimateTform for this tile pair to estimate transform between them
        % add transform types for estimateTforme
        
        registrationSettingsStruct.transformType=transformTypeList;
        if ~isfield(registrationSettingsStruct,'intensityNormType')
            if or(downsampleFactor_fixed~=1,downsampleFactor_move~=1)
                registrationSettingsStruct.intensityNormType='EachImage'; % different resolutions, scale each independently
            else
                registrationSettingsStruct.intensityNormType='SameScaling'; % same resolution. Hopefully same exposure times too
            end
        end
        
        % call esimateTform (my wrapper function for MATLAB's detectSURFFeatures, extractFeatures, matchFeatures, and estimateGeometricTransform2D)
        [tformList,numInliersList,statusFlagList]=estimateTform(imgMove,imgFixed,registrationSettingsStruct); % default is 'rigid'
        
        if ~iscell(tformList)
            tformList={tformList};
        end
        
        fprintf('done.\n')
        % rotation from rigid tform.Rotation
        %temp=padarray(tformList{1}.Rotation,[1 1],'post'); temp(3,3)=3; temp=rotm2eul(temp); rotationFromRigid=180/pi*temp(1); % make a 3x3 and get euler angles
        
        % put translations, rotations (from euler angle in Z), and
        % determinant in lists
        numTforms=length(transformTypeList);
        translationRowColRawList=nan(numTforms,2);
        determinantList=nan(numTforms,1);
        rotationList=nan(numTforms,1);
        for iTform=1:numTforms
            
            transRowColOffsetRaw=-1*[tformList{iTform}.T(3,2), tformList{iTform}.T(3,1)]; % add this to (row,col) of move to get fixed
            eulerAnglesZYX=180/pi()*rotm2eul(tformList{1}.T);
            rotation=eulerAnglesZYX(1);
            determinant=det(tformList{iTform}.T);
                        
            
          % The UmPerPixel of the images that got sent to estimateTransform
          % for registration were downsampled to be the same UmPerPixel.
          % However, the images may have had different height & widths. If
          % so, even if the image centers were identical in stage space,
          % we'd find a translation since the feature coordinates are referenced to
          % the upper left corner, not the centers of the images. So, we
          % need to subtract out this intrinsic offset for images of
          % different sizes.
          intrinsicRowColOffset= 0.5 * [size(imgFixed,1)-size(imgMove,1),size(imgFixed,2)-size(imgMove,2)]; % positive values if imgFixed has more pixels
          transRowColOffset = transRowColOffsetRaw + intrinsicRowColOffset; % add out the effect of different numbers of pixels in the image
          % convert to XY stage offsets
          stageXYOffsetLocal=RowCol2stageXYOffset(transRowColOffset,cameraAngle_Ref,rowDimIsFlippedVsYDim_Ref,umPerPixelOfRegImgs);
          % now account for the fact that these images may not actually have same XY
          % coordinates
          %stageXYOffset=rstageXY_fixed-stageXY_move+stageXYOffsetLocal;
          stageXYOffset=stageXY_move-rstageXY_fixed+stageXYOffsetLocal;
            
            % put in table for output
            iT=iT+1;
            T.tileID_fixed(iT)=tileID_fixed;
            T.tileID_move(iT)=tileID_move;
            T.transformType(iT)=transformTypeList(iTform);
            T.tform(iT)=tformList(iTform); % has to stay in a cell because different object types
            T.numInliers(iT)=numInliersList(iTform);
            T.statusFlag(iT)=statusFlagList(iTform);
            T.rotation(iT)=rotation;
            T.determinant(iT)=determinant;
            T.transRowColRawOffset(iT,1:2)=transRowColOffsetRaw;
            T.intrinsicRowColOffset(iT,1:2)=intrinsicRowColOffset;
            T.transRowColOffset(iT,1:2)=transRowColOffset;
            T.stageXYOffsetLocal(iT,1:2)=stageXYOffsetLocal;
            T.rstageXY_fixed(iT,1:2)=rstageXY_fixed;
            T.stageXY_move(iT,1:2)=stageXY_move;
            T.stageXYOffset(iT,1:2)=stageXYOffset;
            
            % put in array lists for fprintf
            translationRowColRawList(iTform,1:2)= transRowColOffsetRaw;
            rotationList(iTform)=rotation; % could also maybe use rotm2axang(rotm)
            determinantList(iTform)=determinant;
            
        end
        
        
        % print out results of the transforms
        fprintf('                         %22s%22s%22s\n',transformTypeList{1},transformTypeList{2},transformTypeList{3})
        fprintf('    numInliers:         %22i%22i%22i\n',numInliersList(1),numInliersList(2),numInliersList(3))
        fprintf('    statusFlag:         %22i%22i%22i\n',statusFlagList(1),statusFlagList(2),statusFlagList(3))
        fprintf('    translation (pixels):[%9.3f %9.3f] [%9.3f %9.3f] [%9.3f %9.3f]\n',translationRowColRawList(1,1),translationRowColRawList(1,2),translationRowColRawList(2,1),translationRowColRawList(2,2),translationRowColRawList(3,1),translationRowColRawList(3,2))
        fprintf('    rotations (degrees):%22.3f%22.3f%22.3f\n',rotationList(1),rotationList(2),rotationList(3))
        fprintf('    determinants:       %22.3f%22.3f%22.3f\n',determinantList(1),determinantList(2),determinantList(3))
        fprintf('    %s\n',repmat('-',1,86))
        
        
    end % end iCpPairs loop
    
    T.umPerPixelOfRegImgs=repmat(umPerPixelOfRegImgs,numRowsT,1);
    
    % close readers
    reader_fixed.close()
    reader_move.close()
    
    
end