% call_plotScanObject

clear Scan
scanDir='/Volumes/IAND_05/E159part2_scan/scan/scan1_WellA1';
load(fullfile(scanDir,filesep,'ScanObject.mat'))
%%
%plotScanObject(Scan)
% appdesigner(filename) opens the specified .mlapp file in App Designer. If the .mlapp file is not on the MATLAB® path, specify the full path.
appdesigner('IanTutorialApp2.mlapp')

%%
app=IanTutorialApp2(Scan)
