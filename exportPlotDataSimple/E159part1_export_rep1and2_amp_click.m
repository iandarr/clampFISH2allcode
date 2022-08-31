% export tables for data plotted amplification and click plots
clear
TgeneMap=array2table([{'gfp','tmr','alexa','cy'}',{'UBC','ITGA3','FN1','MITF'}'],'VariableNames',{'channelBaseName','gene'});
varsOut={'rep','condID'	'condName'	'AmpRound'	'SpecialClick'	'channel' 'intensities'	'channelBaseName' 'gene' 'expTime'	'normFactor'	'intensitiesNorm'};

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);

%% Import Data
% Replicate 1
extractedDataDirForAllCondRep1=fullfile(parentDir,'paper/extractedData/E159_scan_amp/analysis1_amp',filesep);
TspotsRep1=readtable(fullfile(extractedDataDirForAllCondRep1,filesep,'Tspots4AllNormalized.csv'),'ReadVariableNames', true);
TspotsOutRep1=join(TspotsRep1,TgeneMap);
TspotsOutRep1.rep=repmat(1,height(TspotsOutRep1),1);
TspotsOutRep1=movevars(TspotsOutRep1,'rep','Before',1);
TspotsOutRep1=TspotsOutRep1(:,varsOut);

% Replicate 2
extractedDataDirForAllCondRep2=fullfile(parentDir,'paper/extractedData/E159part1rep2_ampclick',filesep);
TspotsRep2=readtable(fullfile(extractedDataDirForAllCondRep2,filesep,'Tspots4AllNormalized.csv'),'ReadVariableNames', true);
TspotsOutRep2=join(TspotsRep2,TgeneMap);
TspotsOutRep2.rep=repmat(2,height(TspotsOutRep2),1);
TspotsOutRep2=movevars(TspotsOutRep2,'rep','Before',1);
TspotsOutRep2=TspotsOutRep2(:,varsOut);

% combine
TspotsOut=[TspotsOutRep1; TspotsOutRep2];
%% export SourceData

% Export #1: for Figure 1 only (rep1, UBC)
idxFig1only=all([ismember(TspotsOutRep1.condID,[1 2 3 4 5 6]), strcmp(TspotsOutRep1.channelBaseName,'gfp')],2);
TspotsOut_Fig1only=TspotsOutRep1(idxFig1only,:);
delete(fullfile(dirForDataExport,filesep,'SourceData_Fig1.xlsx'))
writetable(TspotsOut_Fig1only,fullfile(dirForDataExport,filesep,'SourceData_Fig1.xlsx'),'Sheet','SourceData_Fig1')

% Export #2:  CSV for Extended Data Fig 4 and Supplementary Fig 3 (both reps, all genes)
delete(fullfile(dirForDataExport,filesep,'SourceData_ExtDataFig4_Amp4genes.csv'))
writetable(TspotsOut,fullfile(dirForDataExport,filesep,'SourceData_ExtDataFig4_Amp4genes.csv')) % one section

% Export #3: Same thing as #2 but XLSX, and split into two side-by-side sections in excel
outdataFile='SourceData_ExtDataFig4_Amp4genes.xlsx';
outdataSheet='Sheet1';
delete(fullfile(dirForDataExport,filesep,outdataFile))
writecell({'Points 1 to 1,000,000 are in columns A to L'},fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet,'Range','A1')
writetable(TspotsOut(1:1000000,:),fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet,'Range','A2')
writecell({'Points 1,000,001 to the end are in columns N to Y'},fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet,'Range','N1')
writetable(TspotsOut(1000001:end,:),fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet,'Range','N2')


%% Export Sample Size N
groupingVars={'rep','condName','gene','expTime'};
DataVars='intensitiesNorm';
ToutN=sortrows(grpstats(TspotsOut(:,[groupingVars,DataVars]),groupingVars,'median',"DataVars",DataVars),{'rep','gene'})
outdataFileForN='N_ExtDataFig4_Amp4genes.xlsx';
writetable(ToutN,fullfile(dirForDataExport,filesep,outdataFileForN))