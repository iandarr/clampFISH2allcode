function [ControlPoints,cameraAngleVsRefAngle_Move,rowDimIsFlippedVsYDim_Move]=registerControlPoints(Scan,fixedRound,moveRound,numControlPointsPerRegistration,registrationSettingsStruct)
% this version of registerControlPoints assumes we can immediately get the
% right control point pairs (no local search upon failure) but can handle
% camera angle and wrong rowDimIsFlippedVsYDim_Move assumptions.

fprintf('%s\n',repmat('-',1,86+4));
fprintf('registerControlPoints: Registering fixedRound=%i, moveRound=%i\n',fixedRound,moveRound)

% get fixedRound's registered stage XY positions and moveRound's (unregistered) stage XY positions
    rstageXY_fixed=Scan.Rounds(fixedRound).fullTiles.rstage.T{:,{'centerX','centerY'}}; % registered stage coordinates of the fixed round
    stageXY_move=   Scan.Rounds(moveRound).fullTiles.stage.T{:,{'centerX','centerY'}}; % (unregistered) stage coordinates of this round
    
    %% First look to see if we controlPointInitialGuesses, and get stageXY_move_WithInitialGuessTform accordingly
    if ~isempty(Scan.Rounds(moveRound).controlPtInitialGuesses)
        tformInitialGuess=Scan.Rounds(moveRound).controlPtInitialGuesses.tform;
        tformInitialGuessFixedRound=Scan.Rounds(moveRound).controlPtInitialGuesses.fixedRound; 
        stageXY_move_MappedOntoInitGuessFixedRoundStageCoord=transformPointsForward(tformInitialGuess,stageXY_move);
        % tformInitialGuessFixedRound is already registered (must be
        % according to logic used to determine registrationRoundPairs
        % matrix) so now that moveRounds's stageXY points are mapped onto
        % it, can use its registration tranformation
        tformRegistrationOfInitialGuessFixedRound=Scan.Rounds(tformInitialGuessFixedRound).stage2rstage.tform;
        stageXY_move_WithInitialGuessTform=transformPointsForward(tformRegistrationOfInitialGuessFixedRound,stageXY_move_MappedOntoInitGuessFixedRoundStageCoord);
    else % no control points, so no tforms needed right now
        stageXY_move_WithInitialGuessTform=stageXY_move;
    end
    
    % find eligible tiles to choose as control points
    WidthHeight_move=Scan.Rounds(moveRound).umPerPixel*[Scan.Rounds(moveRound).fullTileNumCols,Scan.Rounds(moveRound).fullTileNumRows];
    WidthHeight_fixed=Scan.Rounds(fixedRound).umPerPixel*[Scan.Rounds(fixedRound).fullTileNumCols,Scan.Rounds(fixedRound).fullTileNumRows];
    minFractionalAreaOfOverlap=0.5;
    [iTileEligible_move,iTileEligible_fixed,areaOverlapUmSq]=getOverlappingTilePairsFromTwoLists(stageXY_move_WithInitialGuessTform,rstageXY_fixed,WidthHeight_move,WidthHeight_fixed,minFractionalAreaOfOverlap); % only an exact solution if angle of rounds in relation to one another is zero. Or if they are 90 or 270 degrees apart but square tiles. But close enough.
    
    % exclude the first tiles if there are >=3 eligible tiles. This is
    % because the first tiles are likely to have nonrepresentative hysteresis
    if size(iTileEligible_move,1)>=3
        keepIndices=iTileEligible_move~=1 & iTileEligible_fixed~=1;
        iTileEligible_move=iTileEligible_move(keepIndices);
        iTileEligible_fixed=iTileEligible_fixed(keepIndices);
        areaOverlapUmSq=areaOverlapUmSq(keepIndices);
    end
    
    stageXY_move_Eligible_WithInitialGuessTform=stageXY_move_WithInitialGuessTform(iTileEligible_move,:);
    
    %% First select 4 numControlPointsPerRegistration which are spread out, roughly at the corners
    cpTilePairsList=nan(0,2); %[fixed, move]
    % first try to select 4 corners:
    [rectxBound,rectyBound,~]=minboundrect(stageXY_move_Eligible_WithInitialGuessTform(:,1),stageXY_move_Eligible_WithInitialGuessTform(:,2)); % get minimal bounding rectangle
    % this function from John D'Errico (2021). A suite of minimal bounding objects (https://www.mathworks.com/matlabcentral/fileexchange/34767-a-suite-of-minimal-bounding-objects), MATLAB Central File Exchange. Retrieved August 22, 2021.
    rectXY4Bound=[rectxBound(1:4),rectyBound(1:4)];
    [~,indEligible4Corners]=min(pdist2(stageXY_move_Eligible_WithInitialGuessTform,rectXY4Bound),[],1);
    indEligible4Corners=unique(indEligible4Corners)'; % only unique ind. makes it potentially less than 4 corners.
    cpTilePairsList=[iTileEligible_fixed(indEligible4Corners),iTileEligible_move(indEligible4Corners)]; %[fixed, move] convert indEligible --> iTile
    
    %% Next select any additional points which are far away from other points
    if numControlPointsPerRegistration>size(cpTilePairsList,1)
        % then get pairs in addition to the corners
        numAdditionalPoints=numControlPointsPerRegistration - size(cpTilePairsList,1);
        % ind=selectSpreadOutPoints(XYall,numAdditional,indSelected)
        indEligible_New=selectSpreadOutPoints(stageXY_move_Eligible_WithInitialGuessTform,numAdditionalPoints,indEligible4Corners);
        %cpTilePairsList=[iTileEligible_move(indEligible_New),iTileEligible_fixed(indEligible_New)];
        cpTilePairsList=[iTileEligible_fixed(indEligible_New),iTileEligible_move(indEligible_New)]; % 28-nov-2021
        
    elseif numControlPointsPerRegistration<size(cpTilePairsList,1)
        % trim to fewer pairs
        cpTilePairsList=cpTilePairsList(1:numControlPointsPerRegistration,:);
    end
    
    %% call registerControlPointsWCamOrient - assumes camera orientations in Scan.Rounds(moveRound) are correct
    cameraAngleVsRefAngle_Move=Scan.Rounds(moveRound).cameraAngleVsRefAngle; %rowDimIsFlippedVsYDim_Move_STOREOFINPUTVALUE=cameraAngleVsRefAngle_Move;
    rowDimIsFlippedVsYDim_Move=Scan.Rounds(moveRound).rowDimIsFlippedVsYDim;
    
    
    foundCpRegistrationStoppingCondition=false; registrationAttempt=0; ambiguousResultCount=0; % initialize
    while foundCpRegistrationStoppingCondition==false
        registrationAttempt=registrationAttempt+1;
        
        if registrationAttempt>1
            fprintf('%s\n',repmat('-',1,100))
            fprintf('----- cameraOrientation assumption of moveRound=%i detected to be %s. Now attempting another registration ----- \n',moveRound,nextRegistrationAttemptStr)
            fprintf('----- registrationAttempt ----- = %i\n',registrationAttempt)
        end
        
        %% call registerControlPointsWCamOrient
        TCpPairs=registerControlPointsWCamOrient(Scan,fixedRound,moveRound,cpTilePairsList,rstageXY_fixed,stageXY_move,cameraAngleVsRefAngle_Move,rowDimIsFlippedVsYDim_Move,registrationSettingsStruct);
    
        % look at registration details in TCpPairs to determine if it's
        % good
        TCpPairs_Rigid=     TCpPairs(strcmp(TCpPairs.transformType,'rigid'),:); 
        TCpPairs_Similarity=TCpPairs(strcmp(TCpPairs.transformType,'similarity'),:);
        TCpPairs_Affine=    TCpPairs(strcmp(TCpPairs.transformType,'affine'),:);
        assert(height(TCpPairs_Rigid)==size(cpTilePairsList,1))
        
        if median(TCpPairs_Rigid.numInliers)>4
            % then we have a good rigid transform. Check that angle is 0
            rotationOfCpPairsWithInliers=TCpPairs_Rigid.rotation(TCpPairs_Rigid.numInliers>=2); % with 2 or more control points
            rotationEstimate=wrapTo180(middleAngle(rotationOfCpPairsWithInliers));
            rotationEstimateToNearest90=90*round(rotationEstimate/90);
            
            if rotationEstimateToNearest90==0
                % all is good, don't need to attempt another registration
                ControlPoints=TCpPairs_to_ControlPoints(TCpPairs);
                foundCpRegistrationStoppingCondition=true;
            
            elseif rem(rotationEstimateToNearest90,90)<5 % transform is close to a multiple of 90 degrees
                % then we can be pretty sure this should be the camera
                % angle of move round
                cameraAngleVsRefAngle_Move= rotationEstimateToNearest90 + cameraAngleVsRefAngle_Move;
                nextRegistrationAttemptStr=sprintf('wrong cameraAngleVsRefAngle, adjusted to %f',cameraAngleVsRefAngle_Move);
                % now try registration again
            else
                warning('rotation of control points in rigid transform is not close to a multiple of 90 degrees. Not confident in this registration')
                ControlPoints=TCpPairs_to_ControlPoints(TCpPairs);
                foundCpRegistrationStoppingCondition=true; % still output, not sure what to do.
            end
            
        elseif all([median(TCpPairs_Affine.numInliers)>4,...
                    median(TCpPairs_Similarity.numInliers)<median(TCpPairs_Affine.numInliers),...
                    median(TCpPairs_Affine.determinant)<0])
                % this is evidence that rows and columns are flipped
                rowDimIsFlippedVsYDim_Move=~rowDimIsFlippedVsYDim_Move;
                nextRegistrationAttemptStr=sprintf('probably have wrong rowDimIsFlippedVsYDim assumption - this will be flipped to %i (the opposite of previous registration attempt)',rowDimIsFlippedVsYDim_Move);
                % now try registration again
        else
            % ambiguous result
            ambiguousResultCount=ambiguousResultCount+1;
            if ambiguousResultCount==1
                fprintf('Ambiguous registration result. Could be too few matching registration objects, perhaps due to the tile positions not referencing the same area. Or perhaps we have moveRound %i rowDimIsFlippedVsYDim=%i assumption wrong. Will next try reversing this assumption\n',moveRound,rowDimIsFlippedVsYDim_Move);
                
                cameraAngleVsRefAngle_Move_WithFirstAmbiguousResult=cameraAngleVsRefAngle_Move;
                rowDimIsFlippedVsYDim_Move_WithFirstAmbiguousResult=rowDimIsFlippedVsYDim_Move;
                ControlPoints_WithFirstAmbiguousResult=TCpPairs_to_ControlPoints(TCpPairs);
                
                rowDimIsFlippedVsYDim_Move=~rowDimIsFlippedVsYDim_Move;
                nextRegistrationAttemptStr=sprintf('UNKNOWN (registration is not strong)');
            else % second time we got an ambiguous result
                warning('SECOND ambiguous registration result for control points between fixedRound=%i and moveRound=%i. Will return previous ambiguous result')
                cameraAngleVsRefAngle_Move=cameraAngleVsRefAngle_Move_WithFirstAmbiguousResult;
                rowDimIsFlippedVsYDim_Move=rowDimIsFlippedVsYDim_Move_WithFirstAmbiguousResult;
                ControlPoints=ControlPoints_WithFirstAmbiguousResult;
                
                foundCpRegistrationStoppingCondition=true;
            end
        end
    
        if registrationAttempt==4
            error('could not register rounds fixedRound=%i and moveRound=%i',fixedRound,moveRound)
        end
    end
    

end

    
function ControlPoints=TCpPairs_to_ControlPoints(TCpPairs)
        TCpPairs_Rigid=TCpPairs(strcmp(TCpPairs.transformType,'rigid'),:);
        
        numCpPair=height(TCpPairs_Rigid);
        ControlPoints=struct();
        for iCpPair=1:numCpPair
            ControlPoints(iCpPair).fixedRound=TCpPairs_Rigid.fixedRound(iCpPair);
            ControlPoints(iCpPair).tileID_fixed=TCpPairs_Rigid.tileID_fixed(iCpPair);
            ControlPoints(iCpPair).tileID_move=TCpPairs_Rigid.tileID_move(iCpPair);
            
            ControlPoints(iCpPair).fixed_rstage_centerXY=TCpPairs_Rigid.rstageXY_fixed(iCpPair,1:2);
            %ControlPoints(iCpPair).fixed_rstage_centerY=TCpPairs_Rigid.rstageXY_fixed(iCpPair,2);
            
            ControlPoints(iCpPair).move_stage_centerXY=TCpPairs_Rigid.stageXY_move(iCpPair,1:2);
            %ControlPoints(iCpPair).move_stage_centerY=TCpPairs_Rigid.stageXY_move(iCpPair,2);
            
            ControlPoints(iCpPair).move_OffsetXY=TCpPairs_Rigid.stageXYOffset(iCpPair,1:2);
            %ControlPoints(iCpPair).move_OffsetY=TCpPairs_Rigid.stageXYOffset(iCpPair,2);
        end
end