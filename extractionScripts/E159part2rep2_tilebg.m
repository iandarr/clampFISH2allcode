% make single-tile background models for each channel

%% Load Tcond and Tscan
%clear
parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2/paper';
Tcond=readtable(fullfile(parentDir,'experiments/E159part2rep2_scan_Tcond.xlsx'));
Tscan=readtable(fullfile(parentDir,'experiments/E159part2rep2_scan_Tscan.xlsx'));
Tscan.condIDs=cellfun(@(x) str2num(x),Tscan.condIDs,'UniformOutput',false);
cd(parentDir)
scanIDsList=Tscan.scanID;
scanID=3;

% load ScanObject
rowInTscan=find(Tscan.scanID==scanID);
scanDir=Tscan.scanDir{rowInTscan};
load(fullfile(scanDir,filesep,'ScanObject.mat')) % loads Scan struct

backgroundImgDir='/Volumes/IAND_08/rawData/E159part2rep2_scan/scan/scan3_WellA3withStrip/BackgroundTiles';

%
numTilesToSample=100; % use this many tiles to assemble background intensity image
kthSortedToUse=5; % 1 would mean for a given pixel row & column, use the lowest-intensity pixel out of all numTilesToSample tiles. 2 would mean use 2nd-to-lowest intensity.
maxTileIDToConsider=507; % what this 'scan' was composed of
gaussianSigma=50;
rng(0)
sampleTileList=sort(randperm(maxTileIDToConsider,numTilesToSample));

for thisRound=1:Scan.numRounds % each 'round' here is an imaging cycle
    reader = bfGetReader();
    reader = loci.formats.Memoizer(reader,0);
    reader.setId(Scan.Rounds(thisRound).Nd2_filepath);

    channelNames=Scan.Rounds(thisRound).channelNames;
    channelPrefixes=Scan.Rounds(thisRound).channelPrefixes;

    channelExposureTimesMs=Scan.Rounds(thisRound).channelExposureTimesMs;
    channelExposureTimesMsCell=replace(cellstr(num2str(channelExposureTimesMs')),' ','')';
    numChannels=length(channelNames);
    channelOutNames=join([channelPrefixes',repmat({'_'},numChannels,1),channelExposureTimesMsCell',repmat({'ms'},numChannels,1)],'',2)';


    imgCollection=zeros(Scan.Rounds(thisRound).fullTileNumRows,Scan.Rounds(thisRound).fullTileNumCols,numTilesToSample,'uint16');
    for iChannel=1:numChannels
        channelName=channelNames{iChannel};
        %channelName='CY3';
        channelOutName=channelOutNames{iChannel};
        for iTile=1:numTilesToSample
            tileID=sampleTileList(iTile);
            fprintf('%s: getting %s iTile=%3i of %i (tileID=%3i)\n',channelOutName,channelName,iTile,numTilesToSample,tileID)
            img_uint16=getPlaneFromNd2file(reader, tileID, channelName,'closeReaderAfterGettingPlane',false); % added z plane functionality
            imgCollection(:,:,iTile)=img_uint16;
        end
        %imgCollectionSortedAll=sort(imgCollection,3);
        %meansOfSorted=reshape(mean(imgCollectionSortedAll,[1 2]),size(imgCollectionSortedAll,3),1);
        %plot(meansOfSorted)
        imgCollectionSorted=mink(imgCollection,kthSortedToUse,3);
        kthSortedImg=imgCollectionSorted(:,:,kthSortedToUse);

        %     myfittype = fittype('a + b*log(x)',...
        %     'dependent',{'y'},'independent',{'x'},...
        %     'coefficients',{'a','b'});
        %     ranksToFit=1:5;
        %     fitobject=fit(ranksToFit',meansOfSorted(ranksToFit),myfittype)
        %     plot(fitobject)

        % smooth it
        bgImg=imgaussfilt(kthSortedImg,gaussianSigma); % std default is 0.5

        %     % visualize
        %     figure(1); clf
        %     t=tiledlayout(1,2);
        %     ax1=nexttile(t,1); title(sprintf('kthSortedImg, k=%i',kthSortedToUse)); hold on; axis equal;
        %     clims=[100 180];
        %     imagesc(kthSortedImg,clims); colorbar;
        %     ax2=nexttile(t,2); title(sprintf('bgImg (gaussian filtered with sigma=%.2f)',gaussianSigma)); hold on; axis equal;
        %     imagesc(bgImg,clims); colorbar;
        %     linkaxes([ax1,ax2])

        % write result to file
        imwrite(bgImg,fullfile(backgroundImgDir,filesep,[channelOutName,'.tif']))
    end

    reader.close
end