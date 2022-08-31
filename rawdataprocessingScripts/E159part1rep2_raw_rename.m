% E157part2rep2_VsSmFISH rawdataprocessing renaming
%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E159part1rep2_ampclick_Tcond.xlsx'));
parentDirData='/Volumes/IAND_08';
% rename files table
Trename=readtable(fullfile(parentDir,filesep,'paper/experiments/E159part1rep2_ampclick_ImagingDetails.xlsx'));
Trename=Trename(strcmp(Trename.magnification,'100X'),:);
%% rename files
 
condList=Tcond.condID;
assert(isequal(Tcond.condID,[1:height(Tcond)]'))
%magnificationList={'60X','20X','10X'};
magnificationList={'100X'};

cd(parentDirData)

for iMag=1:length(magnificationList)
    magnification=magnificationList{iMag};

    for iCond=1:length(condList)
        condID=condList(iCond);
        tiffDir=Tcond.(['Tiff_',magnification,'_dir']){condID};
        cd(fullfile(parentDirData,tiffDir))
        
        % rename
        
        nRows=height(Trename)
        origList=join([repmat({'img_'},nRows,1),Trename.ElementsInternalChannelName,repmat({'_XY'},nRows,1)],'',2)
        newList=join([Trename.RajlabimagetoolsChannelName,repmat({'0'},nRows,1)],'',2)
%         origList=TrenameBAD.RajlabimagetoolsChannelName;
%         newList=Trename.RajlabimagetoolsChannelName;
%         origList=origList(2:end-1)
%         newList=newList(2:end-1)
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