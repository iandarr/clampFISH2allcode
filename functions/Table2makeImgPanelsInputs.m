function S=Table2makeImgPanelsInputs(T)

numPanels=height(T);
% 
scalingInput=cell(2,numPanels);
%scalingInput(1,:)=repmat({{10,99.9}},1,numPanels);
ch1contrast=[T.ch1ContrastMin,T.ch1ContrastMax];
scalingInput(1,:)=mat2cell(ch1contrast,ones(size(ch1contrast,1),1),2)';

ch2contrast=[T.ch2ContrastMin,T.ch2ContrastMax];
scalingInput(2,:)=mat2cell(ch2contrast,ones(size(ch2contrast,1),1),2)';
% folder and file name
%assert(all(strcmp(Tex.tifFolder,Tex.tifFolder(1))))
%tifFolder=Tex.tifFolder{1};
%imgFiles=cellfun(@(x) fullfile(tifFolder,filesep,[x,'.tif']),[Tex.ch1Name,Tex.ch2Name],'UniformOutput',false)';
imgFiles(1,:)=join([T.tifFolder';T.ch1Name';repmat({'.tif'},1,numPanels)],'',1);
imgFiles(2,:)=join([T.tifFolder';T.ch2Name';repmat({'.tif'},1,numPanels)],'',1);

%segmentationImgs=cellfun(@(x) fullfile([x,'.tif']),Tex.segmentationImg,'UniformOutput',false)';
segmentationImgs=join([T.tifFolder';T.segmentationImg';repmat({'.tif'},1,numPanels)],'',1);
inPanelLabels=T.inPanelLabel';


% boundingbox, layoutMatrix, outColors
boundingbox=cellfun(@(x) str2num(x),T.boundingbox,'UniformOutput',false)';
layoutMatrix=zeros(max(T.panelRow),max(T.panelCol));
ind=sub2ind(size(layoutMatrix),T.panelRow,T.panelCol);
layoutMatrix(ind)=1:numPanels;
outColors=[cellfun(@(x) str2num(x),T.ch1Color,'UniformOutput',false)';...
           cellfun(@(x) str2num(x),T.ch2Color,'UniformOutput',false)'];
if ismember('rotateDegrees',T.Properties.VariableNames)
       rotateDegrees=num2cell(T.rotateDegrees)';
else
    rotateDegrees=num2cell(zeros(1,numPanels));
end

if ismember('zplanes',T.Properties.VariableNames)
    zplanes=(T.zplanes)';
else
    zplanes=repmat({'all'},1,numPanels);
end


% boundingbox may get shifted
if ismember('xShift',T.Properties.VariableNames)
    boundingbox=cellfun(@(x,shift) [x(1)+shift, x(2), x(3), x(4)],boundingbox,num2cell(T.xShift)','UniformOutput',false);
end

if ismember('yShift',T.Properties.VariableNames)
    boundingbox=cellfun(@(x,shift) [x(1), x(2)+shift, x(3), x(4)],boundingbox,num2cell(T.yShift)','UniformOutput',false);
end 

%[imgFiles,scalingInput,segmentationImgs,inPanelLabels,boundingbox,layoutMatrix,outColors]
S.imgFiles=imgFiles;
S.scalingInput=scalingInput;
S.segmentationImgs=segmentationImgs;
S.inPanelLabels=inPanelLabels;
S.boundingbox=boundingbox;
S.layoutMatrix=layoutMatrix;
S.outColors=outColors;
S.rotateDegrees=rotateDegrees;
S.zplanes=zplanes;
end