% E157part2rep2_VsSmFISH rawdataprocessing renaming
%% Table of all conditions

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
cd(parentDir)
% Read in table of conditions (Tcond)
Tcond=readtable(fullfile(parentDir,filesep,'paper/experiments/E157part2rep2_VsSmFISH_Tcond.xlsx'));
parentDirData='/Volumes/IAND_07';
% rename files table
%%
%% rename files
 
condList=Tcond.condID;
%magnificationList={'60X','20X','10X'};
magnificationList={'10X'};

cd(parentDirData)

for iMag=1:length(magnificationList)
    magnification=magnificationList{iMag};

    for iCond=1:length(condList)
        condID=condList(iCond);
        tiffDir=Tcond.(['Tiff_',magnification,'_dir']){condID};
        cd(fullfile(parentDirData,tiffDir)); cd('../'); % channelMapInfoDir
        TchannelMapInfo=readtable('ChannelMapInfo.xlsx');
        cd(fullfile(parentDirData,tiffDir))
        
        % rename
        
        nRows=height(TchannelMapInfo)
        origList=join([repmat({'img_'},nRows,1),TchannelMapInfo.ChannelName_KeyInput,repmat({'_XY'},nRows,1)],'',2)
        newList=join([TchannelMapInfo.Rajlabimagetools_ValueOutput,repmat({'0'},nRows,1)],'',2)
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