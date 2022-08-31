% export simple data for RNAseq experiments
clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/';

dirForDataExport=fullfile(parentDir,filesep,'paper/plotsDataSimple',filesep);
outdataFile='SupFig10_RNAseq.xlsx';
outdataSheet='SupFig10_RNAseq';
delete(fullfile(dirForDataExport,filesep,outdataFile)) % because writing to excel sheet only overwrites in the designated range. If existing data is outside that range it will remain.

% import data
% rep1 = E159part2_scan_VsRNAseq_plot2

Tdataout=table();

for rep=1:2
    if rep==1
        tableDir='paper/extractedData/E159part2_scan/VsRNAseq/';
    elseif rep==2
        tableDir='paper/extractedData/E159part2rep2_scan/VsRNAseq/';
    else
        error('')
    end

    T_RNAseq_naive=readtable(fullfile(parentDir,tableDir,'tpm_RNAseq_naiveSubset.csv'));
    T_RNAseq_resistant=readtable(fullfile(parentDir,tableDir,'tpm_RNAseq_resistantSubset.csv'));
    T_clamp2_naive=readtable(fullfile(parentDir,tableDir,'clamp2_meanCount_naive.csv'));
    T_clamp2_resistant=readtable(fullfile(parentDir,tableDir,'clamp2_meanCount_resistant.csv'));

    T=[T_RNAseq_naive;T_RNAseq_resistant;T_clamp2_naive;T_clamp2_resistant];
    Tdataout_ThisRep=table();
    Tdataout_ThisRep.replicate=[0; 0; rep; rep];
    Tdataout_ThisRep.experimentID={'';'';'E159part2';'E159part2rep2'};
    Tdataout_ThisRep.datatype={'RNAseq';'RNAseq';'clampFISH2';'clampFISH2'};
    Tdataout_ThisRep.units={'tpm';'tpm';'MeanNumberOfSpotsPerCell';'MeanNumberOfSpotsPerCell'};
    Tdataout_ThisRep.celltype={'drug-naive WM989 A6-G3';'vemurafenib-resistant WM989 A6-G3 RC4';'drug-naive WM989 A6-G3';'vemurafenib-resistant WM989 A6-G3 RC4'};
    Tdataout_ThisRep=[Tdataout_ThisRep,T];

    Tdataout=[Tdataout;Tdataout_ThisRep];
end
assert(isequal(Tdataout(1:2,:),Tdataout(5:6,:))) % should be same RNAseq data
Tdataout=Tdataout([1,2,3,4,7,8],:); % remove RNAseq duplicates

writetable(Tdataout,fullfile(dirForDataExport,filesep,outdataFile),'Sheet',outdataSheet)

