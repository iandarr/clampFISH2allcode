% E157part2rep2_VsSmFISH_sensitivity
% clampFISH 2.0 vs. smFISH correlation
%% Import experiment description Tcond
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
%parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable('paper/experiments/E157part2rep2_VsSmFISH_Tcond.xlsx');
if ~all(Tcond.condID==[1:height(Tcond)]')
    error('condID must be equivalent its row number in this code')
end

extractedDataDirForAllCond='paper/extractedData/E157part2rep2_VsSmFISH/';
TcellMap2=readtable(fullfile([extractedDataDirForAllCond,filesep,'TcellMapBetweenMagnifications2.csv']));


condID=1;
spots60X=TcellMap2.numSpots_60X(TcellMap2.condID==condID);
spots20X=TcellMap2.numSpots_20X(TcellMap2.condID==condID);
format compact
% Sensitivity and specificity with DDX58: single threshold

t60X=3; % 'hi' threshold for 60X smFISH
t20X=3; % 'hi' threshold for 20X clampFISH 2.0

TP=sum( spots20X >= t20X & spots60X >= t60X );
TN=sum( spots20X <  t20X & spots60X <  t60X );
FP=sum( spots20X >= t20X & spots60X <  t60X );
FN=sum( spots20X <  t20X & spots60X >= t60X );

sensitivity=100* TP / (TP + FN);
specificity=100* TN / (TN + FP);

numCellsTotal=length(spots60X);
assert((TP + TN + FP + FN )== length(spots20X))

fprintf('\n')
fprintf("sharp 'hi' threshold 60X >=%i and 20X>=%i:\n",t60X,t20X)
fprintf('TotalCells=%i, TP=%i, TN=%i, FP=%i, FN=%i\n',numCellsTotal,TP,TN,FP,FN)
fprintf('specificity=%.1f (%i/%i), sensitivity=%.1f (%i/%i)\n',specificity,TN,(TN + FP), sensitivity,TP,(TP + FN))
