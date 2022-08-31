% E159part2_scan_VsRNAseq_extract
% high-throughput scan, average FISH vs. bulk RNAseq
% for:
%   - drug-naive WM989 A6-G3 (wells A1,A2,A3,B1,B2) 
%   - 1uM vemurafenib-resistant WM989 A6-G3 RC4 (well B3)


parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';
Tcond=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'paper/experiments/E159part2rep2_scan_Tscan.xlsx'));
extractedDataDir=fullfile(parentDir,'paper/extractedData/E159part2rep2_scan',filesep);
outTableDir=fullfile(extractedDataDir,filesep,'VsRNAseq/');

cd(parentDir)

% import Tcell (for drug-naive cells)
Tcell=readtable(fullfile(extractedDataDir,'TcellAll2.csv'));
Tcell.isResistant=logical(Tcell.isResistant);
Tcell_naive=Tcell(~Tcell.isResistant,:);
Tcell_resistant=Tcell(Tcell.isResistant,:);
fprintf('%i drug-naive cells, %i resistant cells (%i total)\n',height(Tcell_naive),height(Tcell_resistant),height(Tcell))
%% find means of FISH data

% use UBC only from round 2
varNames={'R1_CY3_WNT5A','R1_A594_DDX58','R1_CY5_AXL','R2_YFP_UBC','R2_CY3_NGFR','R2_A594_FN1','R2_CY5_EGFR','R3_CY3_ITGA3','R3_A594_MMP1','R3_CY5_MITF'};
temp=split(varNames','_');
geneSymbols=temp(:,3)';

clamp2_meanCount_naive=array2table(mean(Tcell_naive{:,varNames},1),'VariableNames',geneSymbols);
clamp2_meanCount_resistant=array2table(mean(Tcell_resistant{:,varNames},1),'VariableNames',geneSymbols);

writetable(clamp2_meanCount_naive,fullfile(outTableDir,'clamp2_meanCount_naive.csv'))
writetable(clamp2_meanCount_resistant,fullfile(outTableDir,'clamp2_meanCount_resistant.csv'))

%% Load RNAseq data
RNAseqAll=readtable(fullfile(parentDir,filesep,'rawData/RNAseq/wm989bulkSeqYoGo_abv.csv'));
%unique(RNAseqAll.sampleID)
sampleID_naive='YGWM989';
sampleID_resistant='YGRC4PLX1uM';
RNAseq_naive=RNAseqAll(strcmp(RNAseqAll.sampleID,sampleID_naive),:);
RNAseq_resistant=RNAseqAll(strcmp(RNAseqAll.sampleID,sampleID_resistant),:);
%% output tpm summary
RNAseq_naive_10gene=RNAseq_naive(ismember(RNAseq_naive.GeneSymbol,geneSymbols),:);
RNAseq_resistant_10gene=RNAseq_resistant(ismember(RNAseq_resistant.GeneSymbol,geneSymbols),:);

tpm_RNAseq_naiveSubset=array2table((RNAseq_naive_10gene.tpm)','VariableNames',(RNAseq_naive_10gene.GeneSymbol)');
tpm_RNAseq_resistantSubset=array2table((RNAseq_resistant_10gene.tpm)','VariableNames',(RNAseq_resistant_10gene.GeneSymbol)');

% harmonize variable ordering
[~,newOrder] = ismember(geneSymbols,tpm_RNAseq_naiveSubset.Properties.VariableNames);
tpm_RNAseq_naiveSubset=tpm_RNAseq_naiveSubset(:,newOrder);
[~,newOrder] = ismember(geneSymbols,tpm_RNAseq_resistantSubset.Properties.VariableNames);
tpm_RNAseq_resistantSubset=tpm_RNAseq_resistantSubset(:,newOrder);


writetable(tpm_RNAseq_naiveSubset,fullfile(outTableDir,'tpm_RNAseq_naiveSubset.csv'))
writetable(tpm_RNAseq_resistantSubset,fullfile(outTableDir,'tpm_RNAseq_resistantSubset.csv'))

