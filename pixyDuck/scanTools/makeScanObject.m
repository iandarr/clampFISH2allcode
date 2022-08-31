function Scan=makeScanObject(Nd2FilepathList, varargin)
%optional Args:
% ‘referenceRound’,1, (defaults to round 1)
% ‘UmPerPixel’,scalarOrVector, (defaults to finding this in each Nd2 file)
% ‘scanMatrix’,scanMatrix,

%% Input parser
ip=inputParser;

% mandatory inputs
ip.addRequired('Nd2FilepathList',@(x) or(iscell(x),ischar(x)))

% optional inputs, immediately fill in with defaults
ip.addParameter('referenceRound',1,@(x) all([isnumeric(x),isscalar(x)]))
ip.addParameter('controlPointInitialGuessesTable',[],@istable)
ip.addParameter('numControlPointsPerRegistration',9,@(x) and(isnumeric(x),isscalar(x))) % First 4 are the 4 corners (excluding tile 1), additional are the 
ip.addParameter('channelLabels',[],@(x) or(iscell(x),istable(x))) % EITHER: numRoundsx1 cell array of 1xnumChannels cell arrays OR a table with VariableNames round, channel,channelLabels, and any other column with data to be associated with each channel, such as channelPrefixes
ip.addParameter('channelLabelsExactMatchOnly',false,@islogical) % if a channelLabels table is provided: if channelLabelsExactMatchOnly==true, then the name in the channelLabels.channel variable must be an exact match with the channel in the ND2 file. If false, only must start with that channel name. Eg. the properties for Round 1 channels 'YFP' and 'YFP_1' and 'YFP_2' will all get the properties of the row with round=1 and channel='YFP'.
ip.addParameter('skipGetInnerBox',false,@(x) islogical(x)|isnumeric(x))

% optional inputs, determine values later on (leave empty for now)
ip.addParameter('cameraAngleOfReferenceRound',[],@(x) all([isnumeric(x),x>=-180,x<=360]))
ip.addParameter('rowDimIsFlippedVsYDim',[],@(x) and(size(x,2)==1,all(any([x==1,x==0],2))));
ip.addParameter('registrationRoundPairs',[],@(x) and(isnumeric(x),size(x,2)==2))
ip.addParameter('cameraOrientationSettings',struct(),@isstruct)
ip.addParameter('registrationSettings',struct(),@isstruct)

%% Parse inputs
ip.KeepUnmatched=true;
ip.parse(Nd2FilepathList,varargin{:})

% mandatory input
Nd2FilepathList=ip.Results.Nd2FilepathList;

% optional inputs, immediately fill in with defaults
referenceRound=ip.Results.referenceRound;
controlPointInitialGuessesTable=ip.Results.controlPointInitialGuessesTable;
skipGetInnerBox=ip.Results.skipGetInnerBox;

% optional inputs, determine values later on (leave empty for now)
cameraAngleOfReferenceRound=ip.Results.cameraAngleOfReferenceRound;
rowDimIsFlippedVsYDim=ip.Results.rowDimIsFlippedVsYDim;
registrationRoundPairsInput=ip.Results.registrationRoundPairs;
cameraOrientationSettings=ip.Results.cameraOrientationSettings;
numControlPointsPerRegistration=ip.Results.numControlPointsPerRegistration;
registrationSettings=ip.Results.registrationSettings;
channelLabelsInput=ip.Results.channelLabels;
channelLabelsExactMatchOnly=ip.Results.channelLabelsExactMatchOnly;

%% check over inputs
if ischar(Nd2FilepathList)
    Nd2FilepathList=cellstr(Nd2FilepathList); %if single character vector with filename
end

numRounds=numel(Nd2FilepathList);
assert(min(size(Nd2FilepathList))==1) % column or row vector

for thisRound=1:numRounds
    if ~isfile(Nd2FilepathList{thisRound})
        error('could not find Nd2 filepath %s. Is current directory right?',Nd2FilepathList{thisRound})
    end
end

% check that reference round is in range
assert(ismember(referenceRound,1:numRounds))

fprintf('makeScanObject started %s:\n',char(datetime))

%% get a reader list
% % this isn't reliable for some reason. Possibly need to not close the
% reader
% fprintf('  getting bioformats readers for all rounds\n')
% readerList=cell(numRounds,1);
% for thisRound=1:numRounds
%     ND2filepath=Nd2FilepathList{thisRound};
%     fprintf('    round %i: making reader with Nd2 file %s...',thisRound,ND2filepath);
%     isValidReader = @(x) isa(x, 'loci.formats.IFormatReader') && ~isempty(x.getCurrentFile());
%
%     reader = bfGetReader();
%     reader = loci.formats.Memoizer(reader,0);
%     reader.setId(ND2filepath);
%     %reader.close
%     %reader.setId(ND2filepath);
%     fprintf('done.\n');
%     readerList{thisRound}=reader;
%
%     if ~isValidReader(reader)
%         warning('round %i reader is not valid. May need to restart matlab',thisRound)
%     end
%
%     % %reader.close % doesn't work when we close the readers
%     %
% end


%% Initialize Scan object with basic file info
% file info, with all rounds, including the XY positions of each tile

% scan-level stuff
Scan=struct();
Scan.numRounds=numRounds;
Scan.referenceRound=referenceRound;
Scan.cameraAngleOfReferenceRound=cameraAngleOfReferenceRound;
Scan.registrationRoundPairsInput=registrationRoundPairsInput;
Scan.controlPointInitialGuessesTable=controlPointInitialGuessesTable;

%% loop 1 through rounds: add basic info to Scan.Rounds (pre-registration)
for thisRound=1:numRounds
    
    % general Scan info
    Nd2_filepath=Nd2FilepathList{thisRound};
    [~,tmpb,tempc]=fileparts(Nd2_filepath);
    Nd2_filename=[tmpb,tempc];
    Scan.Rounds(thisRound).Nd2_filepath=Nd2_filepath;
    Scan.Rounds(thisRound).Nd2_filename=Nd2_filename;
    
    Scan.Rounds(thisRound).isRef=thisRound==Scan.referenceRound;
    
    
    % old individual functions to get metadata
    % reader=readerList{round}; % doesn't work anymore
    %[Xpos,Ypos]=nd2positions(reader,'closeReaderAfterGettingPlane',false);%get XY positions - STOPPED WORKING WITH READER (unreliable)
    %[Xpos,Ypos]=nd2positions(Nd2_filepath,'closeReaderAfterGettingPlane',true);    % get XY positions
    % Scan.Rounds(thisRound).UmPerPixel=getNd2micronPerPixel(Scan.Rounds(thisRound).Nd2_filepath,'closeReaderAfterGettingPlane',true); % Get UmPerPixel
    % Scan.Rounds(thisRound).channels = readND2Channels2(Nd2_filepath); % get channels
    
    %Scan.Rounds(thisRound).fullTileNumRows=numRows;
    %Scan.Rounds(thisRound).fullTileNumCols=numCols;
    
    % get metadata: fullTileNumRows, fullTileNumCols, umPerPixel, channelNames, channelExposureTimeMs, and XY positions
    metadataStruct=getMetadata(Nd2_filepath,'errorOnAbnormalFlag',true);
    Scan.Rounds(thisRound).numTiles=                metadataStruct.numStacks;
    Scan.Rounds(thisRound).numChannels=             metadataStruct.numChannels;
    Scan.Rounds(thisRound).numZPerTile=             metadataStruct.numZPerStack;
    Scan.Rounds(thisRound).numTimesPerTile=         metadataStruct.numTimesPerStack;
    Scan.Rounds(thisRound).umPerPixel=              metadataStruct.umPerPixel;
    Scan.Rounds(thisRound).fullTileNumRows=         metadataStruct.numImgRows;
    Scan.Rounds(thisRound).fullTileNumCols=         metadataStruct.numImgCols;
    
    Scan.Rounds(thisRound).channelNames=            metadataStruct.channelNames;
    Scan.Rounds(thisRound).channelExposureTimesMs=  metadataStruct.channelExposureTimesMs;
    Scan.Rounds(thisRound).channelLabels=repmat({''},1,Scan.Rounds(thisRound).numChannels);
    
    Scan.Rounds(thisRound).fullTiles.stage.T=table();
    Scan.Rounds(thisRound).fullTiles.stage.T.tileID=[1:Scan.Rounds(thisRound).numTiles]';
    Scan.Rounds(thisRound).fullTiles.stage.T.centerX=metadataStruct.X;
    Scan.Rounds(thisRound).fullTiles.stage.T.centerY=metadataStruct.Y;
    
end % end rounds

%% determine the number of scan resolutions and 
allRoundResRaw=[Scan.Rounds(:).umPerPixel];
allRoundRes=setCloseValuesToMode(allRoundResRaw); % can have issue where resolution is eps(res) away from another resolution value. Just pick most popular one (mode) and use this.
if ~all(allRoundResRaw'==allRoundRes')
    % then detected at least one close but not equal umPerPixel:
    for thisRound=1:numRounds
        Scan.Rounds(thisRound).umPerPixel=allRoundRes(thisRound);
    end
end

Scan.resolutions=sort(unique(allRoundRes,'stable')); % order of size
Scan.numRes=length(Scan.resolutions);

%% add channelLabels if they were provided by user
channelPropTable=formatChannelAssociatedInputs(Scan,channelLabelsInput,channelLabelsExactMatchOnly);
channelProps=channelPropTable.Properties.VariableNames;
for thisRound=1:numRounds
    for iChanProp=1:length(channelProps)
        channelProp=channelProps{iChanProp};
        Scan.Rounds(thisRound).(channelProp)=channelPropTable.(channelProp){thisRound};
    end
end

%% Add Scan.cameraAngleOfReferenceRound and Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim
if  ~isempty(cameraAngleOfReferenceRound)&&~isempty(rowDimIsFlippedVsYDim)
    % then just store user-input values
    Scan.cameraAngleOfReferenceRound=cameraAngleOfReferenceRound;
    Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim=rowDimIsFlippedVsYDim;
elseif xor(isempty(cameraAngleOfReferenceRound),isempty(rowDimIsFlippedVsYDim))
    error('if cameraAngleOfReferenceRound or rowDimIsFlippedVsYDim is user-input, then the other variable must also be input. Usually rowDimIsFlippedVsYDim=true')
else
    % then find camera angle
    %reader=readerList{referenceRound}; % doesn't work sometimes
    Nd2_filepath=Scan.Rounds(referenceRound).Nd2_filepath;
    ImgRowsCols=[Scan.Rounds(referenceRound).fullTileNumRows,Scan.Rounds(referenceRound).fullTileNumCols];
    
    % get reader & file info
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId(Nd2_filepath);
    XYpos=[Scan.Rounds(referenceRound).fullTiles.stage.T.centerX,Scan.Rounds(referenceRound).fullTiles.stage.T.centerY];
    UmPerPixel=Scan.Rounds(referenceRound).umPerPixel;
    
    
    if ~isfield(cameraOrientationSettings,'featuresMatchThreshold')
        cameraOrientationSettings.featuresMatchThreshold=2;
    end
    
    %featureRecognitionSettingsCellArray=namedargs2cell(featureRecognitionSettingsStruct);
    [cameraAngleOfReferenceRound,rowDimIsFlippedVsYDim,cameraAngleList,numInliersList,statusFlagList]=findCameraOrientation(reader,XYpos,UmPerPixel,ImgRowsCols,cameraOrientationSettings); % with reader
    
    fprintf('final output from findCameraOrientation: cameraAngle=%0.03f, rowDimIsFlippedVsYDim=%i\n',cameraAngleOfReferenceRound,rowDimIsFlippedVsYDim)
    reader.close
    
    Scan.cameraAngleOfReferenceRound=cameraAngleOfReferenceRound;
    Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim=rowDimIsFlippedVsYDim;
end

%% Fill in camera orientation details of other non-reference rounds with default values
% will get checked in registerControlPoints function

for thisRound=1:numRounds
    Scan.Rounds(thisRound).cameraAngleVsRefAngle=0;
    Scan.Rounds(thisRound).rowDimIsFlippedVsYDim=Scan.Rounds(Scan.referenceRound).rowDimIsFlippedVsYDim;
    if thisRound==Scan.referenceRound
        Scan.Rounds(thisRound).cameraOrientationVsRefConfirmed=[]; % not relevant % for reference round, we found the angle, and
    else
        Scan.Rounds(thisRound).cameraOrientationVsRefConfirmed=false; % for non-reference round we still need to confirm angle is approximately equal (to nearest 90 degrees), otherwise update stage tiles
    end
end

%% fill in stage tile info (Scan.Rounds(thisRound).fullTiles.rstage .T and .polyvec and tileOrientationConfirmed)
% for now, need to assume cameraAngle and rowDimIsFlippedVsYDim are both
% equal to that of the reference round. Will update this later if it is wrong (but, angle
% in non-reference rounds can only be +/-90 or +/-180 degrees of the
% referenceRound camera angle)

    cameraAngle=            Scan.cameraAngleOfReferenceRound;
    rowDimIsFlippedVsYDim=  Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim;

for thisRound=1:numRounds
    centerXY=Scan.Rounds(thisRound).fullTiles.stage.T{:,{'centerX','centerY'}};
    numRows=Scan.Rounds(thisRound).fullTileNumRows;
    numCols=Scan.Rounds(thisRound).fullTileNumCols;
    UmPerPixel=Scan.Rounds(thisRound).umPerPixel;
    Tvertices=getTileCornerVertices(centerXY,numRows,numCols,UmPerPixel,cameraAngle,rowDimIsFlippedVsYDim);
    Scan.Rounds(thisRound).fullTiles.stage.T=[Scan.Rounds(thisRound).fullTiles.stage.T, Tvertices];
    Scan.Rounds(thisRound).fullTiles.stage.polyvec=polyshapeVector(Tvertices{:,{'X_RowOneColOne','X_RowOneColEnd','X_RowEndColEnd','X_RowEndColOne'}},Tvertices{:,{'Y_RowOneColOne','Y_RowOneColEnd','Y_RowEndColEnd','Y_RowEndColOne'}});
end

%% check over user-input controlPointInitialGuessesTable
checkControlPointInitialGuessesTable(controlPointInitialGuessesTable,Scan.numRounds,Scan.referenceRound)

%%  get registrationRoundPairs matrix
% this is the order in which they'll be registered.
Scan.registrationRoundPairs=chooseRegistrationRoundPairs(Scan.numRounds,Scan.referenceRound,registrationRoundPairsInput,controlPointInitialGuessesTable); % order is important - can't register to a round that isn't already registered itself (unless it's the referenceRound)

%% look at controlPointInitialGuessesTable, add this info to the Scan structure
for thisRound=1:numRounds
    if or(isempty(Scan.controlPointInitialGuessesTable),~istable(controlPointInitialGuessesTable))
        Scan.Rounds(thisRound).controlPtInitialGuesses=[];
    elseif ~any(controlPointInitialGuessesTable.moveRound==thisRound)
        Scan.Rounds(thisRound).controlPtInitialGuesses=[];
    else
        rowsInCpTable=find(controlPointInitialGuessesTable.moveRound==thisRound);
        Scan.Rounds(thisRound).controlPtInitialGuesses.fixedRound=controlPointInitialGuessesTable.fixedRound(rowsInCpTable(1))'; % fixedRound is same for all these rows
        Scan.Rounds(thisRound).controlPtInitialGuesses.fixedX=    controlPointInitialGuessesTable.fixedX(rowsInCpTable)';
        Scan.Rounds(thisRound).controlPtInitialGuesses.fixedY=    controlPointInitialGuessesTable.fixedY(rowsInCpTable)';
        Scan.Rounds(thisRound).controlPtInitialGuesses.moveX=     controlPointInitialGuessesTable.moveX(rowsInCpTable)';
        Scan.Rounds(thisRound).controlPtInitialGuesses.moveY=     controlPointInitialGuessesTable.moveY(rowsInCpTable)';
        
        % get a geometric transform
        xyFixed=[[Scan.Rounds(thisRound).controlPtInitialGuesses.fixedX]',[Scan.Rounds(thisRound).controlPtInitialGuesses.fixedY]'];
        xyMove =[[Scan.Rounds(thisRound).controlPtInitialGuesses.moveX ]',[Scan.Rounds(thisRound).controlPtInitialGuesses.moveY ]'];
        %fitgeotrans(xyMove,xyFixed,transformationType)
        numCps=size(xyFixed,1);
        if numCps==1
            % add on a fake cp that is offset by 1 in X and Y
            xyFixed=[xyFixed;[xyFixed+1 ]];
            xyMove= [xyMove; [xyMove+1  ]];
        end
        if size(xyFixed,1)==2
            tform=fitgeotrans(xyMove , xyFixed, 'nonreflectivesimilarity'); % no reflection
        else
            tform=fitgeotrans(xyMove , xyFixed, 'similarity'); % allows reflection
            if det(tform.T)<0
                warning('the controlPointInitialGuessesTable provided resulted in a coordinate transformation with a reflection')
            end
        end
        Scan.Rounds(thisRound).controlPtInitialGuesses.tform=tform;
    end
end

%% For reference round, set registered stage coordinates (rstage) equal to stage coordinates (stage)
Scan.Rounds(referenceRound).fullTiles.rstage=Scan.Rounds(referenceRound).fullTiles.stage;
Scan.Rounds(referenceRound).stage2rstage.fixedRound=referenceRound; % self
Scan.Rounds(referenceRound).stage2rstage.tform=affine2d; % default identity tform
%% choose control point tiles and register them
% for each non-reference round
% 
% find corresponding control point tiles tbased on controlPtInitialGuessesTform (if they exist) and the stage XY positions

for thisPair=1:size(Scan.registrationRoundPairs,1)
    
    % note thisRound == moveRound
    fixedRound=Scan.registrationRoundPairs(thisPair,1);
    thisRound= Scan.registrationRoundPairs(thisPair,2); % 'move' round (= this round)
    % choose control point tiles. 
    
    [ControlPoints,cameraAngleVsRefAngle_Move,rowDimIsFlippedVsYDim_Move]=registerControlPoints(Scan,fixedRound,thisRound,numControlPointsPerRegistration,registrationSettings);
    % Save control point information
    Scan.Rounds(thisRound).ControlPoints=ControlPoints;
    
    % if registration found camera orientation to be different than saved
    % default, update it
    tileOrientationNeedsUpdating=false;
    if Scan.Rounds(thisRound).cameraAngleVsRefAngle~=cameraAngleVsRefAngle_Move
        fprintf('based on Control Point registration, updating round %i cameraAngleVsRefAngle=%.1f\n',thisRound,cameraAngleVsRefAngle_Move)
        Scan.Rounds(thisRound).cameraAngleVsRefAngle=cameraAngleVsRefAngle_Move;
        tileOrientationNeedsUpdating=true;
    end
    if Scan.Rounds(thisRound).rowDimIsFlippedVsYDim~=rowDimIsFlippedVsYDim_Move
        fprintf('based on Control Point registration, updating round %i rowDimIsFlippedVsYDim=%i\n',thisRound,rowDimIsFlippedVsYDim_Move)
        Scan.Rounds(thisRound).rowDimIsFlippedVsYDim=rowDimIsFlippedVsYDim_Move;
        tileOrientationNeedsUpdating=true;
    end    
    Scan.Rounds(thisRound).cameraOrientationVsRefConfirmed=true;
    
    if tileOrientationNeedsUpdating==true
        % update stage tiles with new cameraOrientation 
        cameraAngle=Scan.cameraAngleOfReferenceRound + Scan.Rounds(thisRound).cameraAngleVsRefAngle;
        Tvertices=getTileCornerVertices(Scan.Rounds(thisRound).fullTiles.stage.T{:,{'centerX','centerY'}},Scan.Rounds(thisRound).fullTileNumRows,Scan.Rounds(thisRound).fullTileNumCols,Scan.Rounds(thisRound).umPerPixel,cameraAngle,Scan.Rounds(thisRound).rowDimIsFlippedVsYDim);
        Scan.Rounds(thisRound).fullTiles.stage.T=[Scan.Rounds(thisRound).fullTiles.stage.T, Tvertices];
        Scan.Rounds(thisRound).fullTiles.stage.polyvec=polyshapeVector(Tvertices{:,{'X_RowOneColOne','X_RowOneColEnd','X_RowEndColEnd','X_RowEndColOne'}},Tvertices{:,{'Y_RowOneColOne','Y_RowOneColEnd','Y_RowEndColEnd','Y_RowEndColOne'}});
    end
    
    % based on control points, fill in
    % Scan.Rounds(tformInitialGuessFixedRound).stage2rstage.tform to go
    % from stageXY --> registered stageXY
    Scan.Rounds(thisRound).stage2rstage.tform=ControlPoints2Tform(Scan.Rounds(thisRound).ControlPoints);
    
    % use tform to generate rstage coordinates
    rstage_centerXY=transformPointsInverse(Scan.Rounds(thisRound).stage2rstage.tform,Scan.Rounds(thisRound).fullTiles.stage.T{:,{'centerX','centerY'}});
    Scan.Rounds(thisRound).fullTiles.rstage.T=array2table([Scan.Rounds(thisRound).fullTiles.stage.T.tileID,rstage_centerXY],'VariableNames',{'tileID','centerX','centerY'});
    cameraAngle=Scan.cameraAngleOfReferenceRound + Scan.Rounds(thisRound).cameraAngleVsRefAngle;
    Tvertices=getTileCornerVertices(rstage_centerXY,Scan.Rounds(thisRound).fullTileNumRows,Scan.Rounds(thisRound).fullTileNumCols,Scan.Rounds(thisRound).umPerPixel,cameraAngle,Scan.Rounds(thisRound).rowDimIsFlippedVsYDim);
    Scan.Rounds(thisRound).fullTiles.rstage.T=[Scan.Rounds(thisRound).fullTiles.rstage.T, Tvertices];
    Scan.Rounds(thisRound).fullTiles.rstage.polyvec=polyshapeVector(Tvertices{:,{'X_RowOneColOne','X_RowOneColEnd','X_RowEndColEnd','X_RowEndColOne'}},Tvertices{:,{'Y_RowOneColOne','Y_RowOneColEnd','Y_RowEndColEnd','Y_RowEndColOne'}}); 
    
%         tileTable.regstagecoord_OffsetX=interp1(...
%         [Scan.Rounds(thisRound).ControlPoints(1:end).stagecoord_CenterX_move],...
%         [Scan.Rounds(thisRound).ControlPoints(1:end).stagecoord_OffsetX],...
%         Scan.Rounds(thisRound).stagecoord_CenterX,'linear','extrap');
%     tileTable.regstagecoord_OffsetY=interp1(...
%         [Scan.Rounds(thisRound).ControlPoints(1:end).stagecoord_CenterY_move],...
%         [Scan.Rounds(thisRound).ControlPoints(1:end).stagecoord_OffsetY],...
%         Scan.Rounds(thisRound).stagecoord_CenterY,'linear','extrap');
end


%% select a stage coordinate to be the origin (at 0,0) for rscan (registered scan) coordinates

Scan.rscanOriginXY=selectScanOrigin(Scan.Rounds(referenceRound).fullTiles.rstage.T,Scan.cameraAngleOfReferenceRound,Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim);

%% prepopulate rscan resolution info. Could have just 1 resolution if all images are same pixel size.
for iRes=1:Scan.numRes
    res=Scan.resolutions(iRes);
    for thisRound=1:numRounds
        if Scan.Rounds(thisRound).umPerPixel==res
            Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).isNativeRes=true;
            Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).nativePxPerRscanPx=1;
        else
            Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).isNativeRes=false;
            Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).nativePxPerRscanPx=res/Scan.Rounds(thisRound).umPerPixel;
        end
    end
end
%% generate rscan.res(1), rscan.res(2),... for fullTile
numRounds=Scan.numRounds;

for iRes=1:Scan.numRes
    for thisRound=1:numRounds
        
        % get centers of rscan fullTiles
        stageCenterXY=Scan.Rounds(thisRound).fullTiles.rstage.T{:,{'centerX','centerY'}};
        umPerPixel_scan=Scan.resolutions(iRes);
        nativePxPerRscanPx=Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).nativePxPerRscanPx;
        numImgRows=Scan.Rounds(thisRound).fullTileNumRows;
        numImgCols=Scan.Rounds(thisRound).fullTileNumCols;
        rowDimIsFlippedVsYDim_refRound=Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim;
        rowDimIsFlippedVsYDim_thisRound=Scan.Rounds(thisRound).rowDimIsFlippedVsYDim;
        cameraAngleVsRefAngle=Scan.Rounds(thisRound).cameraAngleVsRefAngle;
        %                   rstage2rscan(stageXY,rscanOriginXY,cameraAngle,rowDimIsFlippedVsYDim,UmPerPixel,numRowIsEven,numColIsEven,isNativeRes)
        %rscanRowColCenterExact=rstage2rscan(stageXY,      rscanOriginXY,     cameraAngle,                      rowDimIsFlippedVsYDim,                          UmPerPixel)
        rscanRowColCenterExact= rstage2rscan(stageCenterXY,Scan.rscanOriginXY,Scan.cameraAngleOfReferenceRound,rowDimIsFlippedVsYDim_refRound,                  umPerPixel_scan);%,numRowIsEven,numColIsEven,isNativeRes);
        
        % get rscan Table and polyvec
        tileIDs=Scan.Rounds(thisRound).fullTiles.rstage.T.tileID;
        %                                                T=getrscanTable(tileIDs,rscanRowColCentersExact,numImgRows,numImgCols,nativePxPerRscanPx,cameraAngleVsRefAngle)
        T=getrscanTable(tileIDs,rscanRowColCenterExact,numImgRows,numImgCols,nativePxPerRscanPx,cameraAngleVsRefAngle);
        Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).T=T;
        Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(T.TopRow,T.BottomRow,T.LeftCol,T.RightCol);

%         % store box region encompassing all of this round
        Tbox=getOuterBoxFromTilesTable(T);
        Scan.Rounds(thisRound).regions.fullTiles.outerBox.rscan.res(iRes).T=Tbox;
        Scan.Rounds(thisRound).regions.fullTiles.outerBox.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(Tbox.TopRow,Tbox.BottomRow,Tbox.LeftCol,Tbox.RightCol);
        
        % store the largest contained box region associated with this round
        fprintf('getting InnerBox region for round %i \n',thisRound);
        if (skipGetInnerBox==false) && ismember(thisRound,skipGetInnerBox)
            Tbox=getInnerBoxFromPolyvec(Scan.Rounds(thisRound).fullTiles.rscan.res(iRes).polyvec);
            Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res(iRes).T=Tbox;
            Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(Tbox.TopRow,Tbox.BottomRow,Tbox.LeftCol,Tbox.RightCol);
            Scan.Rounds(thisRound).hasInnerBox=true;
        else
            Scan.Rounds(thisRound).hasInnerBox=false;
        end
        %         image coordinates
                                              %T=getFullTileImageTable(tileIDs,numImgRows,numImgCols,cameraAngleVsRefAngle,rowDimIsFlippedVsYDim_refRound,rowDimIsFlippedVsYDim_thisRound)
        %Scan.Rounds(thisRound).fullTiles.image.T=getFullTileImageTable(tileIDs,numImgRows,numImgCols,cameraAngleVsRefAngle,rowDimIsFlippedVsYDim_refRound,rowDimIsFlippedVsYDim_thisRound);
        Scan.Rounds(thisRound).fullTiles.image.T=getUprightFullTileImageTable(tileIDs,numImgRows,numImgCols,cameraAngleVsRefAngle);

    end %end rounds loop
    
end % end res loop
%% for full scan, get intersection of outerBoxes and innerBoxes and union of outerBoxes
for iRes=1:Scan.numRes
    % outerBox
    TouterBoxAllRound=table();
    for thisRound=1:numRounds
        TouterBoxAllRound=[TouterBoxAllRound;Scan.Rounds(thisRound).regions.fullTiles.outerBox.rscan.res(iRes).T];
    end
    % outerBox Intersection
    TboxIntersect=getBoxIntersection(TouterBoxAllRound);
    Scan.regions.fullTiles.outerBoxIntersect.rscan.res(iRes).T=      TboxIntersect;
    Scan.regions.fullTiles.outerBoxIntersect.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(TboxIntersect.TopRow,TboxIntersect.BottomRow,TboxIntersect.LeftCol,TboxIntersect.RightCol);
    % outerBox Union
    TboxUnion=getBoxUnion(TouterBoxAllRound);
    Scan.regions.fullTiles.outerBoxUnion.rscan.res(iRes).T=      TboxUnion;
    Scan.regions.fullTiles.outerBoxUnion.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(TboxUnion.TopRow,TboxUnion.BottomRow,TboxUnion.LeftCol,TboxUnion.RightCol);
    
    % innerBox
    if all(Scan.Rounds(thisRound).hasInnerBox)
        TinnerBoxAllRound=table();
        for thisRound=1:numRounds
            TinnerBoxAllRound=[TinnerBoxAllRound;Scan.Rounds(thisRound).regions.fullTiles.innerBox.rscan.res(iRes).T];
        end
        % get innerBox Intersection
        TboxIntersect=getBoxIntersection(TinnerBoxAllRound);
        Scan.regions.fullTiles.innerBoxIntersect.rscan.res(iRes).T=      TboxIntersect;
        Scan.regions.fullTiles.innerBoxIntersect.rscan.res(iRes).polyvec=RowColStartEnd2polyshape(TboxIntersect.TopRow,TboxIntersect.BottomRow,TboxIntersect.LeftCol,TboxIntersect.RightCol);
    end
end

end