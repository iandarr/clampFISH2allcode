function [fig,scalingMinMaxVector]=makeImgPanels(imgFiles,outColors,varargin)

%% important inputs
ip=inputParser;
% imgFiles
% outColors
ip.addParameter('segmentationImgs',[],@(x) (iscell(x) && (size(x,1)==1)) || ischar(x))

%% inputs related to the whole figure
ip.addParameter('outFile','imagePanel',@ischar);
ip.addParameter('rotateDegrees',[],@iscell);
ip.addParameter('layoutMatrix',[],@(x) isnumeric(x));
ip.addParameter('inPanelLabels',[],@(x) iscell(x) || all(ismember(x,{'none'})));
ip.addParameter('fontSize',8,@(x) validateattributes(x,{'numeric'},{'scalar','>=',0}));
ip.addParameter('segmentationSmoothPixels',10,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
ip.addParameter('segmentationLineWidth',0.5,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
ip.addParameter('segmentationCentermostOnly',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
ip.addParameter('figWidth',7.5*72,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
ip.addParameter('exportFigToFile',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
ip.addParameter('exportRGBimgs',true,@(x)validateattributes(x,{'logical'},{'scalar'}));

%% inputs related to scale bar
ip.addParameter('UmPerPixel',[],@(x) validateattributes(x,{'numeric','scalar','positive'}));
ip.addParameter('scaleBarPosition','LowerLeft',@(x) ismember(x,{'','none','LowerLeft','LowerRight','UpperLeft','UpperRight'}));
ip.addParameter('scaleBarPanels',[],@(x) validateattributes(x,{'numeric','vector'}))


%% inputs to be passed to mergeImgsToRGB
ip.addParameter('zplanes','all',@(x) any([strcmp('all',x),ismatrix(x),iscell(x)]));
ip.addParameter('scalingMinMax',{},@iscell)
ip.addParameter('alpha',[],@iscell)
ip.addParameter('boundingbox',[],@(x) (isnumeric(x) && isequal(size(x),[1 4])) || (iscell(x) && isvector(x)))
ip.addParameter('convertTo','uint8',@(x) all(ismember(lower(x),{'uint16','uint8'})))
ip.parse(varargin{:})

% inputs related to the whole figure, or inherently should be input as cell array of length numImgPanels  
outFile=ip.Results.outFile;
layoutMatrix=ip.Results.layoutMatrix;
inPanelLabels=ip.Results.inPanelLabels;
fontSize=ip.Results.fontSize;
segmentationSmoothPixels=ip.Results.segmentationSmoothPixels;
segmentationLineWidth=ip.Results.segmentationLineWidth;
figWidth=ip.Results.figWidth;
exportFigToFile=ip.Results.exportFigToFile;
exportRGBimgs=ip.Results.exportRGBimgs;

segmentationCentermostOnly=ip.Results.segmentationCentermostOnly;
rotateDegrees=ip.Results.rotateDegrees;

% scale bar related inputs, should be made a cell array:
UmPerPixel=ip.Results.UmPerPixel;
scaleBarPosition=ip.Results.scaleBarPosition;
scaleBarPanels=ip.Results.scaleBarPanels;

% other inputs for this functino, should be made a cell array if not
% provided as such:
segmentationImgs=ip.Results.segmentationImgs;

% inputs for mergeImgsToRGB, should be made a cell array:
zplanes=ip.Results.zplanes;
scalingMinMax=ip.Results.scalingMinMax;
alphaList=ip.Results.alpha;
boundingbox=ip.Results.boundingbox;
convertTo=ip.Results.convertTo;


%% imporant inputs
% imgFiles
imgFiles=checkImgFilesInput(imgFiles);

numImgPanels=size(imgFiles,2);
maxNumImgsPerPanel=size(imgFiles,1);

% related to rotateDegrees
if isempty(rotateDegrees)
    rotateDegrees=num2cell(zeros(1,numImgPanels));
else
    assert(all(cellfun(@(x) isnumeric(x) || isscalar(x),rotateDegrees)))
    assert(length(rotateDegrees)==numImgPanels && isvector(rotateDegrees));
end
for i=1:numImgPanels
    if rem(rotateDegrees{i},90)~=0
        error('rotateDegrees inputs must be multiples of 90')
    end
end


% outColors
if ~iscell(outColors)
    error('outColors must be class cell')
end
if size(outColors,2)==1
    outColors=repmat(outColors,1,numImgPanels);
elseif size(outColors,2)~=numImgPanels
    error('outColors should be size numImgsPerPanelx1 cell array, or a numImgsPerPanel x nPanels cell array')
end

% maskImgs - convert to cell if not already and make sure it's the right length
if ischar(segmentationImgs)
    segmentationImgs={segmentationImgs};
end
if ~isempty(segmentationImgs)
    if length(segmentationImgs)==1
        segmentationImgs=repmat(segmentationImgs,numImgPanels,1);
    elseif length(segmentationImgs)~=numImgPanels
        error('maskImgs should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
    end
end

%% inputs related to the whole figure
% outFile

% layoutMatrix
if isempty(layoutMatrix)
    % default is a row of panels
    layoutMatrix=1:numImgPanels;
else
    % 
    if ~all(ismember(1:numImgPanels,layoutMatrix))
        warning('the provided layoutMatrix does not contain all indices (1 to %i) of the imgFiles provided,',numImgPanels)
    end
end
layoutRows=size(layoutMatrix,1);
layoutCols=size(layoutMatrix,2);

% inPanelLabels
if isempty(inPanelLabels)
    % default to first row, without the .tif
    [~,imgFilenames,~]=fileparts(imgFiles);
    inPanelLabels=imgFilenames(1,:); % first row of images only (assumes DAPI in 2nd row)
elseif all(ismember(inPanelLabels,{'none'}))
    inPanelLabels=repmat({'none'},1,numImgPanels);
else
    assert(length(inPanelLabels)==numImgPanels)
end

% fontSize

%% scale bar related inputs
if ~isempty(UmPerPixel)
    UmPerPixel={UmPerPixel};
end

allImagesAreSameResolution=true;
if isempty(UmPerPixel)
    UmPerPixel=cell(1,numImgPanels);
elseif length(UmPerPixel)==1
    UmPerPixel=repmat(UmPerPixel,numImgPanels,1);
elseif length(UmPerPixel)~=numImgPanels
    error('UmPerPixel should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
else
    if ~all(cellfun(@(x) UmPerPixel{1}==x,UmPerPixel))
        allImagesAreSameResolution=false;
    end
end

scaleBarPosition=cellstr(scaleBarPosition);
if length(scaleBarPosition)==1
    scaleBarPosition=repmat(scaleBarPosition,numImgPanels,1);
elseif length(scaleBarPosition)~=numImgPanels
    error('scaleBarPosition should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
end

if isempty(scaleBarPanels) && allImagesAreSameResolution
    scaleBarPanels=layoutMatrix(end,1); % only lower left panel
elseif isempty(scaleBarPanels) && ~allImagesAreSameResolution
    scaleBarPanels=1:numImgPanels;
else
    assert(all(ismember(scaleBarPanels,1:numImgPanels)))
end

%% inputs to be passed to mergeImgsToRGB
% zplanes
if ischar(zplanes) || isnumeric(zplanes)
    zplanes={zplanes};
end

if length(zplanes)==1
    zplanes=repmat(zplanes,1,numImgPanels);
elseif length(zplanes)~=numImgPanels
    error('zplanes should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
end


% scalingMinMax
if isempty(scalingMinMax)
    scalingMinMax=repmat({{}},1,numImgPanels);
elseif size(scalingMinMax,2)==1
    scalingMinMax=repmat(scalingMinMax,1,numImgPanels);
elseif size(scalingMinMax,2)~=numImgPanels
    error('scalingMinMax should be nImgsPerPanel x 1 cell array, to be replicated for all panels, or nImgsPerPanel x nPanels cell array')
end

% alphaList
if ~isempty(alphaList)
    if length(alphaList)==1
        alphaList=repmat(alphaList,1,numImgPanels);
    elseif length(alphaList)~=numImgPanels
        error('alphaList should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
    end
else
    alphaList=cell(1,numImgPanels);
end

% boundingbox
if isnumeric(boundingbox)
    boundingbox=repmat({boundingbox},1,numImgPanels);
else
    assert(size(boundingbox,2)==numImgPanels)
end

% convertTo
if ischar(convertTo)
    convertTo={convertTo};
end
if isempty(convertTo)
    convertTo=cell(1,numImgPanels);
elseif length(convertTo)==1
    convertTo=repmat(convertTo,1,numImgPanels);
elseif length(convertTo)~=numImgPanels
    error('convertTo should be nPanels x 1 cell array or 1x1 cell array to be replicated for all panels')
end

%% get cropped, contrasted, merged RGB images
imgMergedVector=cell(1,numImgPanels);
scalingMinMaxVector=cell(1,numImgPanels);
for iPanel=1:numImgPanels
    
    imgFilesThisPanel=imgFiles(:,iPanel);
    numImgsThisPanel=length(imgFilesThisPanel);
    outColorsThisPanel=outColors(:,iPanel);
    if length(outColorsThisPanel)>numImgsThisPanel
        outColorsThisPanel=outColorsThisPanel(1:numImgsThisPanel);
    end
    
    [imgMerged,scalingMinMaxList]=mergeImgsToRGB(...
        imgFilesThisPanel,...
        outColorsThisPanel,...
        'zplanes',zplanes{iPanel},...    
        'scalingMinMax',scalingMinMax(:,iPanel),...
        'alpha',alphaList{iPanel},...
        'boundingbox',boundingbox{iPanel},...
        'convertTo',convertTo{iPanel});
    
    if rotateDegrees{iPanel}==0
        imgMergedVector{iPanel}=imgMerged;
    else
        imgMergedVector{iPanel}=rot90(imgMerged,rotateDegrees{iPanel}/90);
    end
    
    scalingMinMaxVector{iPanel}=scalingMinMaxList;
end

%% get segmentation images
if ~isempty(segmentationImgs)

    
    [temp,~,iPanel2iSeg]=unique(cell2table([boundingbox;segmentationImgs]'),'rows');
    segmentationImgsUnique=temp.Var2';
    numSegmentationImgs=length(segmentationImgsUnique);
    % read in segmentation images
    segmentationImgsArray=cell(1,numSegmentationImgs);
    for iSegImg=1:numSegmentationImgs
        iPanel=find(iSegImg==iPanel2iSeg,1);
        thisSegmentationImgName=segmentationImgs{iPanel};
        if ~isfile(thisSegmentationImgName)
            error('could not find segmentation image file %s from directory %s',thisSegmentationImgName,pwd)
        end
        if isempty(boundingbox)
            segmentationImgsArray{iSegImg}=imread(thisSegmentationImgName);
        else
            boundingboxThisPanel=boundingbox{iPanel};
                    PixelRegion={[boundingboxThisPanel(2), boundingboxThisPanel(2)+boundingboxThisPanel(4)-1],...
                        [boundingboxThisPanel(1), boundingboxThisPanel(1)+boundingboxThisPanel(3)-1]};
            segmentationImgsArray{iSegImg}=imread(thisSegmentationImgName,'PixelRegion',PixelRegion);
        end
        
        if rotateDegrees{iPanel}~=0
            segmentationImgsArray{iSegImg}=rot90(segmentationImgsArray{iSegImg},rotateDegrees/90);
        end
    end
    
    % get outlines
    segmentationBoundariesArray=cell(1,numSegmentationImgs);
    for iSegImg=1:numSegmentationImgs 
        iPanel=find(iSegImg==iPanel2iSeg,1);
        % find the segmentation boundaries
        thisSegmentationImg=segmentationImgsArray{iSegImg};
        segIDlist=unique(thisSegmentationImg)';
        segIDlist=segIDlist(segIDlist~=0); % don't include zero
        if segmentationCentermostOnly
            % only want centermost outline
            segCentroidsXY=nan(length(segIDlist),2);
            for iSeg=1:length(segIDlist)
                segID=segIDlist(iSeg);
                %stats=regionprops(thisSegmentationImg==iSegID,'basic');
                [r,c]=find(thisSegmentationImg==segID);
                segCentroidsXY(iSeg,1:2)=[mean(c) mean(r)];
            end
            boundingboxThisPanel=boundingbox{iPanel};
            boundingboxCenterXY=(boundingboxThisPanel(3:4)+1)/2;
            [~,iSegCentermost]=min(vecnorm(segCentroidsXY-boundingboxCenterXY,2,2));
            segIDlist=segIDlist(iSegCentermost);
        end

        Blist={};
        for segID=segIDlist
            B=bwboundaries(thisSegmentationImg==segID,8,'noholes'); % 8 connectivity means diagonals are connected
            % should only be one boundary, but there may be more than one per iSegID
            Blist=[Blist;{B}];
        end
        segmentationBoundariesArray{iSegImg}=Blist;
        segmentationSegIDlistArray{iSegImg}=segIDlist;
    end
end

%% create figure with tiledarray
%fig=figure('visible','off','Units', 'points');
fig=figure('Units', 'points');

figHeight= figWidth * layoutRows/layoutCols;
fig.Position=[1 1 figWidth figHeight]; % [LeftColumn TopRow   WIDTH HEIGHT]
% set(f,'visible','on'); f.Name=fName; 
t=tiledlayout(layoutRows,layoutCols);
%t.PositionConstraint='innerposition';
t.PositionConstraint='outerposition';
t.TileSpacing = 'none';t.Padding = 'tight';
t.TileIndexing='rowmajor';
%t.TileSpacing = 'none';t.Padding = 'tight';


% plot into a figure
for iPanel=1:numImgPanels
    %iAx=iAx+1;
    [tileRow,tileCol]=find(layoutMatrix==iPanel);
    tileNum=sub2ind(size(layoutMatrix'),tileCol,tileRow); % switched since tileindexing is rowmajor
    ax(iPanel)=nexttile(t,tileNum); hold on;
    set(ax(iPanel),'FontSize',fontSize)
    %set(ax(iPanel),'PositionConstraint','innerposition')
    imgMerged=imgMergedVector{iPanel};
    imshow(imgMerged,'Parent',ax(iPanel))
    
    %% add segmentation outlines
    if ~isempty(segmentationImgs)
        iSegImg=iPanel2iSeg(iPanel);
        Blist=segmentationBoundariesArray{iSegImg};
        segIDlist=segmentationSegIDlistArray{iSegImg};
        for iSegID=1:length(Blist)
            B=Blist{iSegID};
            numBoundaries=length(B);
            for iBound=1:numBoundaries
                Bxy=B{iBound};
                if segmentationSmoothPixels>0
                    Bxy=smoothXY(Bxy,segmentationSmoothPixels,5,8);
                end
                Bxy=[Bxy;Bxy(1,:)]; % include first point again
                plot(Bxy(:,2),Bxy(:,1),'w', 'LineWidth', segmentationLineWidth)
            end
        end
    end
    
    %% add inPanelLabels
    inPanelLabel=inPanelLabels{iPanel};
    if ~strcmp(inPanelLabel,'none')
        %xlm=get(ax(iPanel),'XLim');
        %ylm=get(ax(iPanel),'YLim');
        %xText=xlm(1)+diff(xlm)*0.05;
        %yText=ylm(1)+diff(ylm)*0.05;
        %text(xText,yText,inPanelLabel,'FontSize',fontSize,'Color','w','Interpreter','none','FontWeight','bold')
        text(0.05,0.92,inPanelLabel,'Units','normalized','FontSize',fontSize,'Color','w','Interpreter','none','FontWeight','bold')
    end
    
    
    %% export Raster images to subfolder
    if exportRGBimgs
        strRowCol=['r',num2str(tileRow),'c', num2str(tileCol)];
        
        % output folder
        [outDirParent,outPanelName,~]=fileparts(outFile);
        imFileNameOut=[strRowCol,' ',inPanelLabel];
        imFileNameOut=replace(imFileNameOut,{newline,char(13)},' ');
        outSubfolderForImgs=fullfile(outDirParent,filesep,outPanelName,filesep);
        if iPanel==1
            if isfolder(outSubfolderForImgs)
                rmdir(outSubfolderForImgs,'s')
            end
            mkdir(outSubfolderForImgs)
        end
        imFilePath=fullfile(outSubfolderForImgs,[imFileNameOut,'.tif']);
        imwrite(imgMerged,imFilePath)
    end

end % end iPanel loop

linkaxes(ax)
%if nargout==0
set(fig,'visible','on')
%end

% %% export file
if exportFigToFile
[~,~,outExt]=fileparts(outFile);
if isempty(outExt)
    outFile=[outFile,'.pdf'];
end
exportgraphics(fig, outFile, 'ContentType', 'vector','Resolution',300);
end
end


function imgFiles=checkImgFilesInput(imgFiles)
if ischar(imgFiles)
    imgFiles={imgFiles};
end

assert(ismatrix(imgFiles))
for iPanel=1:size(imgFiles,2)
    numImgsThisPanel=sum(cellfun(@(x) ~isempty(x),imgFiles(:,iPanel)));
    for iImg=1:numImgsThisPanel
        imgFile=imgFiles{iImg,iPanel};
        if ~isfile(imgFile)
            error('could not find file %s from directory %s',imgFile,pwd)
        end
    end
end

end


function xySmooth=smoothXY(xy,d,upsampleFactor,minInterpPoints)
    if nargin<=1
        d=40;
    end
    if nargin<=2
        upsampleFactor=5;
    end
    if nargin<=3
       minInterpPoints=8;
    end
        % using polyshape
%         warning('off','MATLAB:polyshape:repairedBySimplify')
%         poly=polyshape(xy(:,1)',xy(:,2)');
%         warning('on','MATLAB:polyshape:repairedBySimplify')
%         xy2=poly.Vertices;
%         segDists=vecnorm(xy2([2:end,1],:)-xy2(1:end,:),2,2);
%         PM=poly.perimeter;

        % without using polyshape
        xy2=xy;
        segDists=vecnorm(xy2([2:end,1],:)-xy2(1:end,:),2,2);
        PM=sum(segDists);
    
    % 
    tOrig=cumsum(segDists([end,1:end-1]))'-segDists(end);
    
    % one vertex every d pixel distance, create 
    %d=40;
    numVertex3=ceil(PM/d);
    if numVertex3<minInterpPoints
        numVertex3=minInterpPoints;
    end
    t3=linspace(tOrig(1),tOrig(end),numVertex3);
    xy3=interp1(tOrig,xy2,t3,'linear');
    
    % create xy points based on spline
    %upsampleFactor=5;
    
    t4=linspace(t3(1),t3(end),length(t3)*upsampleFactor);
    xy4=interp1(t3,xy3,t4,'pchip');
    
    
%     f=figure; ax=gca; hold on; axis equal;
%     %plot(ax,xy(:,1),xy(:,2),'.','MarkerEdgeColor','r'); 
%     plot(ax,xy2(:,1),xy2(:,2),'.','MarkerEdgeColor','b');
%     plot(ax,xy3(:,1),xy3(:,2),'o','MarkerEdgeColor','c');
%     plot(ax,xy4(:,1),xy4(:,2),'.-','LineWidth',0.5)
%     
    xySmooth=xy4;
end
