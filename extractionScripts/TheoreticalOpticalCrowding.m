%  calculate theoretical maximum number of spots per area

parentDir='/Users/iandardani/Dropbox (RajLab)/Shared_IanD/clamp2';
outTableDir=fullfile(parentDir,filesep,'/paper/extractedData/TheoreticalOpticalCrowding');

displayImages=false;
displayImagesMinimumInputSpots=50;

thresh=0.04;

ImgSideLength_px=30; % pixels

spotNumVect=[repmat(1,1,20),repmat(2,1,10),repmat(5,1,4),repmat(10,1,2),15:5:100,125,150,175,200,225,250,300,350,400,450,500,600,700,800,900,1000,2000];

replicates=20; % per element of spot number vector
%warning('changed replicates to lower')

% calculate spot sigma
% Gaussian approximations of fluorescence microscope point-spread function models
% Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin
% https://opg-optica-org.proxy.library.upenn.edu/ao/fulltext.cfm?uri=ao-46-10-1819&id=130945#t7
%
% lambda
% filter set name	    related dyes / FPs	emission filter	    lambda (nm) at midpoint
% DAPI	                DAPI	Chroma      D460/50m	        460
% YFP	                Atto 488, GFP, YFP	Chroma HQ535/30m	535
% CY3	                Cy3, TMR	        Chroma HQ567/15m	567
% A594	                Alexa 594	        Omega 630DF30	    630
% CY5	                Cy5, Atto 647N	    Chroma HQ667/30m	667
%

%lambdaEmission=[535];
lambdaEmission=[667];
%lambdaEmission=[630];
%lambdaEmission=[567];

pixelSizeImagePlane=6500; % nm 6500 is with no binning. 13000 is with 2x2 binning.

% % Objective = 100X/1.45NA
% magnification=100;
% NA=1.45;
% n=1.518; % Refractive index. Nikon Type N immersion oil MXA22165

% % Objective = 60X/1.4NA
% magnification=60;
% NA=1.4;
% n=1.518; % Refractive index. Nikon Type N immersion oil MXA22165

% Objective = 20X/0.75NA
magnification=20;
NA=0.75;
n=1.0003; % Refractive index of Air.
dataOutTableName=sprintf('TheorDetectEfficiencyWith_mag=%ix_NA=%.2f_lambda=%i_thresh=%.3f_px=%.1f_simulateddata.csv',magnification,NA,lambdaEmission,thresh,pixelSizeImagePlane/1000);
summaryOutTableName=sprintf('TheorDetectEfficiencyWith_mag=%ix_NA=%.2f_lambda=%i_thresh=%.3f_px=%.1f_summary.csv',magnification,NA,lambdaEmission,thresh,pixelSizeImagePlane/1000);


pixelSize=pixelSizeImagePlane/magnification; %nm/px

% paraxial optics
theoreticalGaussianFitStd_Paraxial_nm=0.21*lambdaEmission/NA; % not used...
% nonparaxial optics
alpha=asin(NA/n);
k=(2*pi()./lambdaEmission); % wavenumber [1/nm]

theoreticalGaussianFitStd_NonParaxial_nm=...
    (1./(n*k))*...
    ((4 - 7*cos(alpha)^(3/2) + 3*cos(alpha)^(7/2))/...
    (7*(1 - cos(alpha)^(3/2)))       )^(-1/2);

assert(abs(theoreticalGaussianFitStd_NonParaxial_nm/theoreticalGaussianFitStd_Paraxial_nm -1)<0.1) % otherwise likely made mistake

sigmaTheoretical_px=theoreticalGaussianFitStd_NonParaxial_nm/pixelSize; % pixels

%warning('increasing spot size artificially')
%sigmaTheoretical_px=1.3*theoreticalGaussianFitStd_NonParaxial_nm/pixelSize; % pixels
% generate spots within a area

ImgSideLength_nm=ImgSideLength_px*pixelSize;

AreaOfImg_umSq=(ImgSideLength_nm/1000)^2;

Tio=table();
numSimsTot=length(spotNumVect)*replicates;
Tio.numSpotsIn=reshape(repmat(spotNumVect,replicates,1),numSimsTot,1);
Tio.numSpotsOut=nan(height(Tio),1);

rng(0);
sigmaForSpotfinding=0.4;
%aTrousMinThreshFactor=1.5;
fineMeshPerPx=5;
%gaussKernelWidth_fine=11; % pick odd number
gaussKernelWidth_fine=round(6*fineMeshPerPx*sigmaTheoretical_px/2)*2+1; % choose odd number encompassing 6 standard deviations (3 SD on either side)

assert(rem(gaussKernelWidth_fine-1,2)==0)

sigma_fine=sigmaTheoretical_px * fineMeshPerPx;

h0 = fspecial('gaussian',gaussKernelWidth_fine,sigma_fine);
h0xfine=repmat(-(gaussKernelWidth_fine-1)/2:(gaussKernelWidth_fine-1)/2,gaussKernelWidth_fine,1);
h0yfine=repmat((-(gaussKernelWidth_fine-1)/2:(gaussKernelWidth_fine-1)/2)',1,gaussKernelWidth_fine);
h0linear=reshape(h0,[gaussKernelWidth_fine^2,1]);
% rng(0)
fineArrayNumel1D=ImgSideLength_px*fineMeshPerPx;

for iSim=1:numSimsTot
    fineArray=zeros(fineArrayNumel1D);

    numSpotsIn=Tio.numSpotsIn(iSim);
    X=ImgSideLength_px*rand(numSpotsIn,1) + 0.5;
    Y=ImgSideLength_px*rand(numSpotsIn,1) + 0.5;
    X_fine=round(fineMeshPerPx*X);
    Y_fine=round(fineMeshPerPx*Y);

    for iSpot=1:numSpotsIn
        x_fine=X_fine(iSpot);
        y_fine=Y_fine(iSpot);
        hxfine=reshape(x_fine+h0xfine,[gaussKernelWidth_fine^2, 1]);
        hyfine=reshape(y_fine+h0yfine,[gaussKernelWidth_fine^2, 1]);
        idxInbounds=all([[hxfine,hyfine]>=1,[hxfine,hyfine]<=fineArrayNumel1D],2);
        hxfineIn=hxfine(idxInbounds);
        hyfineIn=hyfine(idxInbounds);
        h0linearIn=h0linear(idxInbounds);
        indhfine=sub2ind([fineArrayNumel1D,fineArrayNumel1D],hxfineIn,hyfineIn);
        fineArray(indhfine)=fineArray(indhfine)+h0linearIn;
    end

    fineArrayReshape=reshape(fineArray,[fineMeshPerPx,fineArrayNumel1D^2/fineMeshPerPx]);
    intermediateArrayLinear=sum(fineArrayReshape,1);
    intermediateArrayReshape1=reshape(intermediateArrayLinear,[ImgSideLength_px,fineMeshPerPx*ImgSideLength_px])';
    intermediateArrayReshape2=reshape(intermediateArrayReshape1,[fineMeshPerPx,fineArrayNumel1D^2/(fineMeshPerPx^2)]);
    intermediateArrayReshape3=sum(intermediateArrayReshape2,1);
    imgArray=reshape(intermediateArrayReshape3,ImgSideLength_px,ImgSideLength_px)';

        [Xout,Yout,intensitiesOut]=d2utils.findSpotsaTrous2(imgArray,'sigma',sigmaForSpotfinding,'minRegionalMax',thresh);
        numSpotsOut=length(intensitiesOut);
        Tio.numSpotsOut(iSim)=numSpotsOut;
            if displayImages & (numSpotsIn>=displayImagesMinimumInputSpots)
        f=figure(1);
        t=tiledlayout(f,1,2);
        ax=nexttile(t,1); hold off; 
        imagesc(fineArray); colorbar; axis equal; hold on; title(sprintf('fine array:%3i spots input',numSpotsIn))
        ax.XLim=[0.5 fineArrayNumel1D+0.5]; ax.YLim=[0.5 fineArrayNumel1D+0.5];
        ax=nexttile(t,2); hold off; 
        imagesc(imgArray);colorbar; axis equal; hold on; title(sprintf('simulated image: %3i spot calls (%.1f%%)',numSpotsOut,100*numSpotsOut/numSpotsIn))
        ax.XLim=[0.5 ImgSideLength_px+0.5]; ax.YLim=[0.5 ImgSideLength_px+0.5];
        plot(Yout,Xout,'o','MarkerSize',20,'Color','red')
        %linkaxes(findobj(t.Children,'Type','Axes'))
            end

end
Tio.pctDetected=100*Tio.numSpotsOut./Tio.numSpotsIn;
Tio.numSpotsInPerSqUm=Tio.numSpotsIn/AreaOfImg_umSq;
Tio.numSpotsOutPerSqUm=Tio.numSpotsOut/AreaOfImg_umSq;
TioSummary=grpstats(Tio,'numSpotsIn')
pctDetectedQuery=[100 99 98 95 90 85 80 75 70 65 60 50 45 40 35 30 25 20 15 10 5 2 1]';
%%
[~,indPctUnique]=unique(TioSummary.mean_pctDetected); % need this if have too few replicates (or multiple lines at 100% detection) and some of the numbers end up identical
%indPctUnique=1:height(TioSummary);
x=TioSummary.mean_pctDetected(indPctUnique);
ySpotsOut=TioSummary.mean_numSpotsOutPerSqUm(indPctUnique);
ySpotsIn=TioSummary.mean_numSpotsInPerSqUm(indPctUnique);
vq_spotsOut = interp1(x,ySpotsOut,pctDetectedQuery);
vq_spotsIn = interp1(x,ySpotsIn,pctDetectedQuery);
%vq_spotsOut = interp1(TioSummary.mean_pctDetected,TioSummary.mean_numSpotsOutPerSqUm,pctDetectedQuery)
%vq_spotsIn = interp1(TioSummary.mean_pctDetected,TioSummary.mean_numSpotsInPerSqUm,pctDetectedQuery)
Tinterp=table();
Tinterp.pctDetected=pctDetectedQuery;
Tinterp.spotsActualPerUm2=vq_spotsIn;
Tinterp.spotsDetectedPerUm2=vq_spotsOut;
Tinterp.pctRecalculated=100*Tinterp.spotsDetectedPerUm2./Tinterp.spotsActualPerUm2;
Tinterp
% save to extracted data dir 

writetable(Tio,fullfile(outTableDir,filesep,dataOutTableName))
writetable(Tinterp,fullfile(outTableDir,filesep,summaryOutTableName))