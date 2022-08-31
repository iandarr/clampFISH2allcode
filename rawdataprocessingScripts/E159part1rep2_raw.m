% E159part1rep2 rawdataprocessing
%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E159part1rep2_ampclick_Tcond.xlsx'));
Timaging=readtable(fullfile(parentDir,filesep,'paper/experiments/E159part1rep2_ampclick_ImagingDetails.xlsx'));
parentDirData='/Volumes/IAND_08';
% rename files table

origStr=Timaging.ElementsInternalChannelName;
nRows=height(Timaging);
origStr=join([repmat({'img_'},nRows,1),origStr,repmat({'_XY'},nRows,1)],'',2);
newStr=join([Timaging.RajlabimagetoolsChannelName,repmat({'0'},nRows,1)],'',2);

% origStr2=repmat({'XY'},nRows,1);
% newStr2=repmat({'0'},nRows,1);

Trename=table();
Trename.imagingSession=Timaging.imagingSession;
Trename.magnification=Timaging.magnification;
Trename.origStr=origStr;
Trename.newStr=newStr;
%% rename files
 
condList=Tcond.condID;
magnificationList={'60X'};

cd(parentDirData)

for iMag=1:length(magnificationList)
    magnification=magnificationList{iMag};

    for iCond=1:length(condList)
        condID=condList(iCond);
        cd(Tcond.(['Tiff_',magnification,'_dir']){condID})
        
        % rename
        %tiffFiles=listFilesWithExtension('.tif');
        
        origList=Trename.origStr(strcmp(Trename.magnification,magnification));
        newList=Trename.newStr(strcmp(Trename.magnification,magnification));
        nConversions=length(origList);
        assert(length(origList)==length(newList));
        
        for iConvertPair=1:nConversions
            
            origStr=origList{iConvertPair};
            newStr=newList{iConvertPair};
            
            filesWithOrigStrStruct=dir([origStr,'*']);
            for iFile=1:length(filesWithOrigStrStruct)
                filenameOrig=filesWithOrigStrStruct(iFile).name;
                filenameNew=replace(filenameOrig,origStr,newStr);
                movefile(filenameOrig,filenameNew)
                
            end
        end
        cd(parentDirData)
    end
end