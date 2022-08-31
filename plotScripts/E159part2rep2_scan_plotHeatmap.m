% E159part2rep2_scan_plotHeatmap
% high-throughput scan


parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tscan.xlsx'));
extractedDataDir=fullfile(parentDir,'paper/extractedData/E159part2rep2_scan',filesep);
plotDir=fullfile(parentDir,filesep,'paper/plots/E159part2rep2_scan',filesep);

cd(parentDir)

% import Tcell (for drug-naive cells)
Tcell=readtable(fullfile(extractedDataDir,'TcellAll2.csv'));
Tcell.isResistant=logical(Tcell.isResistant);
fprintf('subsetting only the drug-naive cells (%i out of %i original cells)\n',sum(~Tcell.isResistant),height(Tcell))
Tcell=Tcell(~Tcell.isResistant,:);


% Export simple data used for plot
outdataFile='SupFig16_clustering_rep2.xlsx';
outdataSheet='SupFig16_clustering_rep2';

%% plot heatmap with matlab's clustergram (too slow to do all cells this way)

% SUBSET: choose a way to subset the data to run this. For reference, 8299 cells takes ~85 seconds
%cellSubsetMethod='random';
% cellSubsetMethod='subregion';
cellSubsetMethod='OneOrMoreHighMarkers'; % Have 42,802 cells in this 'OneOrMoreHighMarkers' population. Clustering took over a day to finish this.

useSinglePrecision=false;

useOneUBCround=2; % 'median' use median UBC, 0 = use all of them, 1 to use R1 UBC only, 2 to use R2 UBC only, 3 to use R3 UBC only

% for cellSubsetMethod='random';
%numRandomCells=6653;
numRandomCells=42802;

% for cellSubsetMethod='subregion';
subregionSubset=[44,45,46,47,48];

% for cellSubsetMethod='OneOrMoreHighMarkers';
fishChannels=   {'R1_YFP_UBC','R1_CY3_WNT5A','R1_A594_DDX58','R1_CY5_AXL','R2_YFP_UBC','R2_CY3_NGFR','R2_A594_FN1','R2_CY5_EGFR','R3_YFP_UBC','R3_CY3_ITGA3','R3_A594_MMP1','R3_CY5_MITF'};
countThresholds=[nan           15             10             25           nan           30           100            5           nan           50            40             nan];
assert(isequal(size(fishChannels),size(countThresholds)))



genesAlongColumn=true; % I like true better

fontSize=20;

RowPDist='correlation';
ColumnPDist='correlation';

LogTrans=false;
OptimalLeafOrder=true;


if strcmp(cellSubsetMethod,'subregion')
    idxRowSubset=ismember(Tcell.subregionInd,subregionSubset);
elseif strcmp(cellSubsetMethod,'random')
    rng(0);
    idxRowSubset=false(height(Tcell),1);
    idxRowSubset(randperm(height(Tcell),numRandomCells))=true;
elseif strcmp(cellSubsetMethod,'OneOrMoreHighMarkers')
    idxFishWithThresh=~isnan(countThresholds);
    countThresholdsNotNan=countThresholds(idxFishWithThresh);
    channelsFishWithThresh=fishChannels(idxFishWithThresh);
    idxRowSubset=any(Tcell{:,channelsFishWithThresh}>=countThresholdsNotNan,2);
end

%
a=Tcell{:,channelsFishWithThresh}>=countThresholdsNotNan;
TthreshCounts=array2table(countThresholdsNotNan,'VariableNames',channelsFishWithThresh)
TnumAboveThresh=array2table(sum(a,1),'VariableNames',channelsFishWithThresh)
TpctAboveThresh=array2table(100*sum(a,1)/size(a,1),'VariableNames',channelsFishWithThresh)

% Export simple data used for plot
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
TcellForOut=Tcell;
TcellForOut.IncludedInHeatmap=idxRowSubset;
TcellForOut=movevars(TcellForOut,'IncludedInHeatmap','Before',1);
writetable(TcellForOut,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)
%%

%rowSubset=(5000:10000)'
%rowAXLhigh=find(Tcell.R1_CY5_AXL>prctile(Tcell.R1_CY5_AXL,99));
%rowSubset=[rowSubset;rowAXLhigh];

Tcell2=Tcell(idxRowSubset,:);
numCells=height(Tcell2);

%dataVars=Tcell2.Properties.VariableNames(~ismember(Tcell2.Properties.VariableNames,{'cellID'    'scanID'    'scanWell'    'subregionInd'    'subregionRow'    'subregionCol'    'segID'}));
fishChannels={'R1_YFP_UBC','R1_CY3_WNT5A','R1_A594_DDX58','R1_CY5_AXL','R2_YFP_UBC','R2_CY3_NGFR','R2_A594_FN1','R2_CY5_EGFR','R3_YFP_UBC','R3_CY3_ITGA3','R3_A594_MMP1','R3_CY5_MITF'};

if ischar(useOneUBCround) && strcmp(useOneUBCround,'median')
    error('median still needs work')
    UBC=median(Tcell2{:,fishChannelsAllUBC},2);
elseif useOneUBCround>0
    isUBC=cellfun(@(x) contains(x,'UBC'),fishChannels);
    fishChannelsAllUBC=fishChannels(isUBC);
    fishChannelsAllUBC_Round=cellfun(@(x) str2num(x(2)),fishChannelsAllUBC);
    fishChannelUBC=fishChannelsAllUBC{fishChannelsAllUBC_Round==useOneUBCround};
    UBC=Tcell2{:,fishChannelUBC};
    fishChannelsWithoutUBC=fishChannels(~isUBC);
    countsData=Tcell2{:,fishChannelsWithoutUBC};
    countsData=[countsData,UBC];
    fishChannels=[fishChannelsWithoutUBC,{fishChannelUBC}];
end

countsData=Tcell2{:,fishChannels};
temp=split(fishChannels','_'); temp(:,2)={''};
isUBC=cellfun(@(x) contains(x,'UBC'),fishChannels); temp(~isUBC,[1,2])={''};
fishLabels=join(temp,' ');


cellIDs=Tcell2.cellID;

pctKeep=100*numCells/height(Tcell);
fprintf('making clustergram with %i cells (%.2f%% of all %i  cells) starting %s...\n',numCells,pctKeep,height(Tcell),char(datetime));
%cgObj=clustergram(1+countsData,'LogTrans',true,'ColumnLabels',dataLabels,'RowLabels',cellIDs,'Cluster',1,'Standardize','Row'); % requires bioinformatics toolbox
%cgObj=clustergram(countsData,'ColumnLabels',dataLabels,'RowLabels',cellIDs,'Cluster','all'); % requires bioinformatics toolbox

tic1=tic;

%%

if genesAlongColumn==false
    Standardize='Column';
    ColumnLabels=fishLabels;
    RowLabels=cellIDs;
else
    countsData=countsData'; % transpose data
    Standardize='Row';
    ColumnLabels=cellIDs;
    RowLabels=fishLabels;    
end

if LogTrans==true
    countsData=countsData+1;
end

if useSinglePrecision
    countsData=single(countsData);
end

DisplayRatio=[0.1 0.05]; % Ratio of space that row and column dendrograms occupy. when genesAlongColumn==true then its [TopSideDendro leftsideDendro]

cgObj=clustergram(countsData,...
    'ColumnLabels',ColumnLabels,...
    'RowLabels',RowLabels,...
    'Cluster','all',...
    'Standardize',Standardize,...
    'LogTrans',LogTrans,...
    'RowPDist',RowPDist,...
    'ColumnPDist',ColumnPDist,...
    'Colormap',redbluecmap,...
    'OptimalLeafOrder',OptimalLeafOrder,...
    'DisplayRatio',DisplayRatio); % requires bioinformatics toolbox
elapsedTime=toc(tic1);
fprintf('done.\n Clustering and producing clustergram took %.0f seconds\n',elapsedTime)

titleStr=sprintf('%i cells (subset=%s), Standardize=%s, RowPDist=%s, ColPDist=%s, LogTrans=%i, optOrder=%i',numCells,cellSubsetMethod,Standardize,RowPDist,ColumnPDist,LogTrans,OptimalLeafOrder);
addTitle(cgObj,titleStr);
%ax=gca;
%set(ax,'FontSize',fontSize)


get(cgObj) % show clustergram object


fprintf('Saving clustergram at %s...\n', char(datetime))
save(fullfile(plotDir,filesep,sprintf('cgObj_%s',titleStr)),'cgObj')
fprintf('done at %s\n',char(datetime));
%
fprintf('Rendering with plot...\n')

figNum=9;
fName=titleStr;
if useOneUBCround==0
    fName=[fName,' with all UBC'];
end
delete(figure(figNum))
f = figure(figNum); clf(figNum); f.Name=fName; f.Position=[500 500 1500 400]; f.Color=[1 1 1]; set(f,'visible','on'); %ax=gca;
%set(f, 'Color','w')
%set(ax, 'Color','w')

plot(cgObj,f)
fprintf('done.\n')
%title(f,titleStr,'FontSize',fontSize,'FontWeight','bold')
%
exportgraphics(f, fullfile(plotDir,[f.Name,'.eps']), 'ContentType', 'vector');
%% Time to cluster with clustergram
% rowLenList=[1000 2000 4000 8000 16000];
% clusterTime=nan(length(rowLenList),1);
% for i=1:length(rowLenList)
%     idxRowSubset=1:rowLenList(i);
%     Tcell2=Tcell(idxRowSubset,:);
%     cellIDs=Tcell2.cellID;
%     countsData=Tcell2{:,dataVars};
%     fprintf('making clustergram starting %s for %i rows ...\n',char(datetime),rowLenList(i));
%     tic1=tic;
%     cgObj=clustergram(countsData,'ColumnLabels',fishLabels,'RowLabels',cellIDs,'Cluster','all','Standardize',Standardize,'LogTrans',LogTrans,'RowPDist',RowPDist,'ColumnPDist',ColumnPDist,'Colormap',redbluecmap,'OptimalLeafOrder',OptimalLeafOrder); % requires bioinformatics toolbox
%     elapsedTime=toc(tic1);
%     clusterTime(i)=elapsedTime
% end
% clusterTime
% addTitle(cgObj,sprintf('%i cells, Standardize=%s, RowPDist=%s, ColPDist=%s, LogTrans=%i, optOrder=%i',numCells,Standardize,RowPDist,ColumnPDist,LogTrans,OptimalLeafOrder));

%% Clustering with binarized low/high expression states



