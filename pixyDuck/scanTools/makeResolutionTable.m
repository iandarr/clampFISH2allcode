function [resTable,iResStitch2iResScan]=makeResolutionTable(Scan,resolutionsInput)
resTable=table();
row=0; % initialize
if ~ischar(resolutionsInput)
    erorr("input resolutions can only accept default 'native' or 'all' at this time")
end


if strcmp(resolutionsInput,'native')
    for thisRound=1:Scan.numRounds
        for iResScan=1:Scan.numRes
            if Scan.Rounds(thisRound).fullTiles.rscan.res(iResScan).isNativeRes
                row=row+1;
                isNativeRes=Scan.Rounds(thisRound).fullTiles.rscan.res(iResScan).isNativeRes;
                resTable{row,{'thisRound','iResScan','res','isNativeRes'}}=[thisRound, iResScan, Scan.resolutions(iResScan), isNativeRes];
                assert(Scan.resolutions(iResScan)==Scan.Rounds(thisRound).umPerPixel)
            end
        end
    end
    
elseif strcmp(resolutionsInput,'all')
    for thisRound=1:Scan.numRounds
        for iResScan=1:Scan.numRes
            row=row+1;
                isNativeRes=Scan.Rounds(thisRound).fullTiles.rscan.res(iResScan).isNativeRes;
                resTable{row,{'thisRound','iResScan','res','isNativeRes'}}=[thisRound, iResScan, Scan.resolutions(iResScan), isNativeRes];
        end
    end
end

% resTable.iResScan(6)=8;
% resTable.iResScan(1)=4;
% resTable.iResScan(4)=4;
% resTable.iResScan(7)=2;

%% add on iResStich - one for each unique iResScan
[iResScanUnique,~,idx]=unique(resTable.iResScan);
iResStitchOrdered=[1:length(iResScanUnique)]';
resTable.iResStitch=iResStitchOrdered(idx);
resTable=movevars(resTable,'iResStitch','After','thisRound');

%% iResStitch2iResScan map container
iResStitch2iResScan=containers.Map(iResStitchOrdered,iResScanUnique);
end