%%
% Get a scan object in the variable Scan

%%
figure(26); clf
hold on; axis equal; 
ax1=gca;


roundsToPlot=[1,2,3,4];
roundColors={'r','g','b','c','k'};
coordSys='rscan';
iRes=1;

if contains(coordSys,'scan')
    set(ax1,'YDir','reverse')
end 

for iRound=1:length(roundsToPlot)
    thisRound=roundsToPlot(iRound);
    
    tileIDs=Scan.Rounds(thisRound).fullTiles.stage.T.tileID;
    
    labels=arrayfun(@num2str,tileIDs,'UniformOutput',0);
    roundColor=roundColors{iRound};
    
    if contains(coordSys,'stage')
    % plot polyvec
        plot(Scan.Rounds(thisRound).fullTiles.(coordSys).polyvec);
        
    XY_RowOneColOne=[Scan.Rounds(thisRound).fullTiles.(coordSys).T.X_RowOneColOne , Scan.Rounds(thisRound).fullTiles.(coordSys).T.Y_RowOneColOne];
    XY_Centers=[Scan.Rounds(thisRound).fullTiles.(coordSys).T.centerX,Scan.Rounds(thisRound).fullTiles.(coordSys).T.centerY];
    plot(XY_RowOneColOne(:,1),XY_RowOneColOne(:,2),'o')
    text(XY_Centers(:,1),XY_Centers(:,2),labels,'Color',roundColor)
    text(0.05,0.05+0.05*(iRound-1),sprintf('thisRound=%i',thisRound),'Color',roundColor)

    elseif contains(coordSys,'scan')
        % plot polyvec
            plot(Scan.Rounds(thisRound).fullTiles.(coordSys).res(iRes).polyvec);
            
            ColRow_Centers=[Scan.Rounds(thisRound).fullTiles.(coordSys).res(iRes).T.CenterCol,Scan.Rounds(thisRound).fullTiles.(coordSys).res(iRes).T.CenterRow];
         text(ColRow_Centers(:,1),ColRow_Centers(:,2),labels,'Color',roundColor)
    end
end
title(sprintf('coordinateSystem=%s',coordSys))

% outerBoxIntersect
%plot(Scan.regions.fullTiles.outerBoxIntersect.rscan.res(iRes).polyvec)
%%
%cameraAngle=Scan.cameraAngleOfReferenceRound;
%cameraAngle=269.99
%cameraAngle=270
%rscanOriginXY=selectScanOrigin(Scan.Rounds(referenceRound).fullTiles.(coordSys).T,cameraAngle,Scan.Rounds(referenceRound).rowDimIsFlippedVsYDim);
%plot(rscanOriginXY(1),rscanOriginXY(2),'+','Color','r','MarkerSize',20)
%%
thisRound=3;
plot(Scan.Rounds(thisRound).fullTiles.(coordSys).polyvec)

%%
%polyvec1=Scan.Rounds(1).fullTiles.stage.polyvec;
%tform=Scan.Rounds(1).
%polyvec2=transformPointsForward(tform,polyvec1)

% move_OffsetXY =
%   3×2 single matrix
%   457.4098 -589.1564
%   456.7100 -193.6202
%   457.7751 -193.6663

%% plot control points
figure(15); clf
hold on; axis equal; 
roundsToPlot=[1];
for iRound=1:length(roundsToPlot)
    thisRound=roundsToPlot(iRound);
    
    CP=Scan.Rounds(thisRound).ControlPoints;
    numCpPairs=length(CP);
    
    
    
    labels=arrayfun(@num2str,[1:numCpPairs]','UniformOutput',0);
    XY=reshape([CP(:).move_stage_centerXY],2,numCpPairs)';
    plot(XY(:,1),XY(:,2),'.','Color','b')
    text(XY(:,1),XY(:,2),labels,'Color','b')
    
    offsetXY=reshape([CP(:).move_OffsetXY],2,numCpPairs)';
    XY2=XY+offsetXY;
    plot(XY2(:,1),XY2(:,2),'.','Color','g')
    text(XY2(:,1),XY2(:,2),labels,'Color','g')
    
    % confirm tform works
    
    tform=Scan.Rounds(thisRound).stage2rstage.tform;
    XY2_tform=transformPointsForward(tform,XY)
    XY2
    tformError=XY2_tform-XY2
end


%% confirm tform works









