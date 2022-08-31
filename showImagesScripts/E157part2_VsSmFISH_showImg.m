% E157_smFISHCorrelation exampleImages
%%
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
%parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
parentDataDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';


% rep2
Tcond=readtable('paper/experiments/E157part2rep2_VsSmFISH_Tcond.xlsx');
extractedDataDirForAllCond='paper/extractedData/E157part2rep2_VsSmFISH/';
exampleImagesDir='paper/exampleImages/E157part2rep2_clampVsSmFISH';

condID=1; % condID=1 is DDX58
parentDataDir='/Volumes/IAND_07/';
arrObjNums=[23 6]; % example image 60X Array and Object numbers


geneName=Tcond.ReadoutSetName{Tcond.condID==condID};
outFileName60X=['smFISH_60X_',geneName,'_Cy3'];
if strcmp(geneName,'DDX58')
    dye='Alexa594';
else
    dye='Atto647N';
end
outFileName20X=['clampFISH_20X_',geneName,'_',dye];
outFileName10X=['clampFISH_10X_',geneName,'_',dye];


if ~all(Tcond.condID==[1:height(Tcond)]')
    error('condID must be equivalent its row number in this code')
end

TcellMap2=readtable(fullfile([extractedDataDirForAllCond,filesep,'TcellMapBetweenMagnifications2.csv']));


TcellMap2(TcellMap2.condID==condID,{'condID','arrayNum_60X','objNum_60X','numSpots_60X','numSpots_20X','numSpots_10X','arrayNum_20X','objNum_20X','arrayNum_10X','objNum_10X'})

rowToInspect=find(all([TcellMap2.arrayNum_60X==arrObjNums(1),TcellMap2.objNum_60X==arrObjNums(2),TcellMap2.condID==condID],2))
assert(length(rowToInspect)==1)

TcellMap2(rowToInspect,{'condID','arrayNum_60X','objNum_60X','numSpots_60X','numSpots_20X','numSpots_10X','arrayNum_20X','objNum_20X','arrayNum_10X','objNum_10X'})

condIDcheck=TcellMap2.condID(rowToInspect);
assert(condID==condIDcheck)

mag60X='60X';
mag20X='20X';
mag10X='10X';

cd(parentDataDir)



dir60X=Tcond.(['Tiff_',mag60X,'_dir']){condID};
dir20X=Tcond.(['Tiff_',mag20X,'_dir']){condID};
dir10X=Tcond.(['Tiff_',mag10X,'_dir']){condID};

channel60X=TcellMap2.(['channel_',mag60X]){rowToInspect};
channel20X=TcellMap2.(['channel_',mag20X]){rowToInspect};
channel10X=TcellMap2.(['channel_',mag10X]){rowToInspect};

arrayNum60X=TcellMap2.(['arrayNum_',mag60X])(rowToInspect);
arrayNum20X=TcellMap2.(['arrayNum_',mag20X])(rowToInspect);
arrayNum10X=TcellMap2.(['arrayNum_',mag10X])(rowToInspect);

objNum60X=TcellMap2.(['objNum_',mag60X])(rowToInspect);
objNum20X=TcellMap2.(['objNum_',mag20X])(rowToInspect);
objNum10X=TcellMap2.(['objNum_',mag10X])(rowToInspect);

fileSuffix60X=TcellMap2.(['imageFileSuffix_',mag60X])(rowToInspect);
fileSuffix20X=TcellMap2.(['imageFileSuffix_',mag20X])(rowToInspect);
fileSuffix10X=TcellMap2.(['imageFileSuffix_',mag10X])(rowToInspect);


threshold60X=Tcond.(['spots_',mag60X,'_threshold'])(condID);
threshold20X=Tcond.(['spots_',mag20X,'_threshold'])(condID);
threshold10X=Tcond.(['spots_',mag10X,'_threshold'])(condID);


if ~isfolder(exampleImagesDir)
    mkdir(exampleImagesDir)
end


%imgDir=Tcond.(['Tiff_60X_dir']){condID};
%imgDirdir=Tcond.(['Tiff_20X_dir']){condID};

% for 60X image
%zplanesInput='all' % all planes (always max merge)
zplanesInput=3:12; % alwways a max merge of these
[imgFISHSource,croppedMask,spotRowColZ60X]=getRajlabimagetoolsCroppedCellImage(dir60X,arrayNum60X,objNum60X,channel60X,threshold60X,zplanesInput,true);
[imgDAPISource,~,~]=getRajlabimagetoolsCroppedCellImage(dir60X,arrayNum60X,objNum60X,'dapi',threshold60X,zplanesInput,false);
imgMerged60X=mergeImgsToRGB({imgFISHSource;imgDAPISource},{'w';'b'},'scalingMinMax',{{10 99.9};{10 99.5}},'alpha',[0.7;1]);


fh=figure(1); clf; hold on;
set(fh,'Position',[1500 200 600 600],'Units','points')
imshow(imgMerged60X); ax=gca; 
plot(ax,spotRowColZ60X(:,2),spotRowColZ60X(:,1),'o','MarkerSize',6,'Color',[0 1 0],'LineWidth',0.2)
v=visboundaries(croppedMask,'Color',[1 1 1],'LineWidth',0.15);
imwrite(imgMerged60X,fullfile(parentDir,filesep,exampleImagesDir,filesep,[outFileName60X,'.tif']))%without spots
exportgraphics(gca,fullfile(parentDir,filesep,exampleImagesDir,filesep, [outFileName60X,'.eps']), 'ContentType', 'vector'); % with spots



% for 20X image
zplanesInput='all';
[imgFISHSource,croppedMask,spotRowColZ]=getRajlabimagetoolsCroppedCellImage(dir20X,arrayNum20X,objNum20X,channel20X,threshold20X,zplanesInput,true);
[imgDAPISource,~,~]=                    getRajlabimagetoolsCroppedCellImage(dir20X,arrayNum20X,objNum20X,'dapi',threshold20X,zplanesInput,false);
imgMerged20X=mergeImgsToRGB({imgFISHSource;imgDAPISource},{'w';'b'},'scalingMinMax',{{10 99.95};{10 99.5}},'alpha',[0.7;1]);

fh=figure(2); clf; hold on;
set(fh,'Position',[1500 600 200 200],'Units','points')
imshow(imgMerged20X); ax=gca;
plot(ax,spotRowColZ(:,2),spotRowColZ(:,1),'o','MarkerSize',5,'Color',[1 0 0],'LineWidth',0.1)
v=visboundaries(croppedMask,'Color',[1 1 1],'LineWidth',0.05);
imwrite(imgMerged20X,fullfile(parentDir,filesep,exampleImagesDir,filesep,[outFileName20X,'.tif'])) 
exportgraphics(gca,fullfile(parentDir,filesep,exampleImagesDir,filesep, [outFileName20X,'.eps']), 'ContentType', 'vector'); % with spots

% for 10X image
zplanesInput='all';
[imgFISHSource,croppedMask,spotRowColZ]=getRajlabimagetoolsCroppedCellImage(dir10X,arrayNum10X,objNum10X,channel10X,threshold10X,zplanesInput,true);
[imgDAPISource,~,~]=                    getRajlabimagetoolsCroppedCellImage(dir10X,arrayNum10X,objNum10X,'dapi',threshold10X,zplanesInput,false);
imgMerged10X=mergeImgsToRGB({imgFISHSource;imgDAPISource},{'w';'b'},'scalingMinMax',{{15 99.95};{10 95.5}},'alpha',[0.7;1]);

fh=figure(3); clf; hold on;
set(fh,'Position',[1500 1200 100 100],'Units','points')
imshow(imgMerged10X); ax=gca;
plot(ax,spotRowColZ(:,2),spotRowColZ(:,1),'o','MarkerSize',2,'Color',[1 0 0],'LineWidth',0.05)
v=visboundaries(croppedMask,'Color',[1 1 1],'LineWidth',0.01);
imwrite(imgMerged10X,fullfile(parentDir,filesep,exampleImagesDir,filesep,[outFileName10X,'.tif'])) 
exportgraphics(gca,fullfile(parentDir,filesep,exampleImagesDir,filesep, [outFileName10X,'.eps']), 'ContentType', 'vector'); % with spots