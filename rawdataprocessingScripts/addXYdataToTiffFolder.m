% add XYdata.csv file to Tiff folders
cd(parentDir)

% % E159 part 1
%Tcond=readtable('experiments/E159part1_ampclick_Tcond.xlsx');
%magnificationList={'10X','20X'};

% % E159 part 2
Tcond=readtable('experiments/E159part2_scan_Tcond.xlsx');
magnificationList={'20X'};

condList=Tcond.condID;
%magnificationList={'60X','10X','20X'};
%magnificationList={'20X'};
magnificationList={'10X'};


for iMag=1:length(magnificationList)
    magnification=magnificationList{iMag}
    
    for i=1:length(condList)
        
        condID=condList(i)
        rowInTcond=find(Tcond.condID==condID);
        assert(isscalar(rowInTcond));
        
        Nd2_dir=Tcond.(['Nd2_',magnification,'_dir']){rowInTcond};
        
        % get ND2 XY positions
        cd(Nd2_dir)
        nd2fileList=listFilesWithExtension('.nd2');
        
        Nd2_filenametext=Tcond.(['Nd2_',magnification,'_filenametext']){rowInTcond}
        Nd2_filepath=nd2fileList(contains(nd2fileList,Nd2_filenametext));
        assert(length(Nd2_filepath)==1)
        Nd2_filepath=Nd2_filepath{1};
        
        [Xposition,Yposition]=nd2positions(Nd2_filepath);
        cd(parentDir)
        
        % output a csv table in tiffDir
        tiffDir=Tcond.(['Tiff_',magnification,'_dir']){rowInTcond};
        Zstack=[1:length(Xposition)]';
        TXYdata=array2table([Zstack,Xposition,Yposition],'VariableNames',{'Zstack','Xposition','Yposition'});
        writetable(TXYdata,fullfile(tiffDir,filesep,'XYdata.csv'))
        cd(parentDir)
    end
end