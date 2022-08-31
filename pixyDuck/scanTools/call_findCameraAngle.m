% call findCameraAngle

%% get images
Nd2_filepath='/Volumes/IAND_05/E159part2_scan/rawdata/s2_6well_cycle1/20210602_213621_459__WellA1_ChannelDAPI,YFP,CY3,A594,CY5_Seq0001.nd2';

%tileID_img1_list=[1:5:39,40:5:78,79:5:117,118:5:156,157]
tileID_img1_list=[1,2,3,39,40,77,78,79,117,118]
%tileID_img1_list=94
cameraAngle_list=[];
statusFlag_list=[];
for iTileID=1:length(tileID_img1_list)
    tileID_img1=tileID_img1_list(iTileID)
%tileID_img1=195 %78 my program obviously fails here on edge when StageYDimensionIsNotFlipped=true, cameraAngle =  -179.7177;
%tileID_img1=10 %7 a problem since inliers are these weirdos in the middle


% statusFlag
%   0: No error.
%   1: matchedPoints1 and matchedPoints2 do not contain enough points.
%   2: Not enough inliers have been found
% when statusFlag is 1 or 2, cameraAngle=nan and imgcoord_tform=nan

cameraAngle,statusFlag
cameraAngle_list(iTileID,1)=cameraAngle;
statusFlag_list(iTileID,1)=statusFlag;
end
[tileID_img1_list',cameraAngle_list],statusFlag_list
%% playing around
scanDir=Tscan.scanDir{Tscan.scanID==1};
Scan=load(fullfile(scanDir,filesep,'ScanInfo.mat'),'Scan');
Scan=Scan.Scan;
stagecoordXYpos=Scan.Rounds(1).stagecoordXYpos;
angles=atan2(diff(stagecoordXYpos(:,2)),diff(stagecoordXYpos(:,1)))*180/pi()
%angles=atan([diff(stagecoordXYpos(:,2))./diff(stagecoordXYpos(:,1))])*180/pi();
%angles(diff(stagecoordXYpos(:,2))<0)=angles(diff(stagecoordXYpos(:,2))<0)+180; % since arctan only goes -90 to 90 degrees
angles
%%
figure(1); hold off;
plot(stagecoordXYpos(:,1),stagecoordXYpos(:,2),'.')
hold on; 
set(gca, 'YDir','reverse')
set(gca, 'XDir','reverse')
xlabel('stage X')
ylabel('stage Y')
offset=100;
grid on;
labels=arrayfun(@num2str,[1:size(stagecoordXYpos,1)]','UniformOutput',0);
text(offset+stagecoordXYpos(:,1),offset+stagecoordXYpos(:,2),labels,'FontSize',6)