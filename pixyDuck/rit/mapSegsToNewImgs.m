function [out_dImageFileSuffix,out_dBoundingBoxes,out_dMasksCropped,out_dMasks]=mapSegsToNewImgs(sImageFileSuffixes,sBoundingBoxes,sMasksCropped,sXimgs,sYimgs,dXimgs,dYimgs,sImgSize,dImgSize,sPixelSizeMicrons,dPixelSizeMicrons,XoffsetsMicron,YoffsetsMicron,RowIsDim,ColIsDim)
% mapSegsToNewImgs
% s=source
% d=destination
%
% matlab image coordinates are (row,col)
%
%% Done in E157_analysis2
% Tcell1=extractCells('returnMasks');
numObj=length(sImageFileSuffixes);
assert(size(sBoundingBoxes,1)==numObj)
assert(size(sBoundingBoxes,2)==4)
assert(size(sMasksCropped,1)==numObj)

if or(isnan(XoffsetsMicron),isnan(YoffsetsMicron))
   error('inputs XoffsetsMicron and YoffsetsMicron should be numbers, not nan')
end
if length(XoffsetsMicron)==1
    XoffsetsMicron=repmat(XoffsetsMicron,numObj,1);
end

if length(YoffsetsMicron)==1
    YoffsetsMicron=repmat(YoffsetsMicron,numObj,1);
end

assert(length(XoffsetsMicron)==numObj)
assert(length(YoffsetsMicron)==numObj)




% check whether RowIsDim and ColIsDim are good
signCheck=startsWith({RowIsDim,ColIsDim},{'+','-'});
if ~all(signCheck==1)
    error('make sure RowIsDim and ColIsDim both start with either + or -. For example RowIsDim=+y and ColIsDim=-x')
end
xCheck=endsWith({RowIsDim,ColIsDim},'x');
yCheck=endsWith({RowIsDim,ColIsDim},'y');
xyChecktotal=sum([xCheck;yCheck],1);
if ~all([all(xyChecktotal==1),sum(xCheck)==1,sum(yCheck==1)])
    error('For RowIsDim and ColIsDim character array inputs, one of them should end in x and the other should end in y')
end
%%
out_dImageFileSuffix=nan(numObj,1);
out_dBoundingBoxes=nan(numObj,4);
out_dMasksCropped=cell(numObj,1);
out_dMasks=cell(numObj,1);


s_RowColOfImgCenter=round(sImgSize/2);
d_RowColOfImgCenter=round(dImgSize/2);
for iObj=1:numObj

    % get source info for this source object
    sImageFileSuffix=sImageFileSuffixes(iObj);
    sXimg=sXimgs(sImageFileSuffix);
    sYimg=sYimgs(sImageFileSuffix);
    
    sBbox=sBoundingBoxes(iObj,:);% [colStart, rowStart, colWidth, rowHeight]
    sMaskCropped=sMasksCropped{iObj,:};
    
    XoffsetMicron=XoffsetsMicron(iObj);
    YoffsetMicron=YoffsetsMicron(iObj);
    
    
    % calculate SOURCE X and Y in microns for the center of the bounding
    % box. Call this the s_XcenterPtMicrons and s_YcenterPtMicrons
    s_objCenterPtRowCol=[sBbox(2)+sBbox(4)/2, sBbox(1)+sBbox(3)/2];
    
    s_objCenterPtRowCol_OffsetFromImgCenter = s_objCenterPtRowCol - s_RowColOfImgCenter;
    
    % centerpoint of source bounding box in XY dimensions
    s_XObjcenterPtMicron=nan;
    s_YObjcenterPtMicron=nan;
    switch RowIsDim
        case '+y' % should be this for E157
            s_YObjcenterPtMicron= sYimg + s_objCenterPtRowCol_OffsetFromImgCenter(1) * sPixelSizeMicrons;
        case '-y'
            s_YObjcenterPtMicron= sYimg - s_objCenterPtRowCol_OffsetFromImgCenter(1) * sPixelSizeMicrons;
        otherwise
            error('cannot handle RowIsDim=%s',RowIsDim)
    end
    
    switch ColIsDim
        case '+x'
            s_XObjcenterPtMicron= sXimg + s_objCenterPtRowCol_OffsetFromImgCenter(2) * sPixelSizeMicrons;
        case '-x' % should be this for E157
            s_XObjcenterPtMicron= sXimg - s_objCenterPtRowCol_OffsetFromImgCenter(2) * sPixelSizeMicrons;
        otherwise
            error('cannot handle ColIsDim=%s',ColIsDim)
    end
    % apply X and Y offsets (in microns) to get destination centerPt
    d_XObjcenterPtMicron=s_XObjcenterPtMicron + XoffsetMicron;
    d_YObjcenterPtMicron=s_YObjcenterPtMicron + YoffsetMicron;
    
    % find the destination image whose center is closest to the center of
    % the bounding box
    distsToDestImgs= [[dXimgs - d_XObjcenterPtMicron].^2 + [dYimgs - d_YObjcenterPtMicron].^2].^0.5;
    [~,indDestImg]=min(distsToDestImgs);
    d_XImg=dXimgs(indDestImg);
    d_YImg=dYimgs(indDestImg);
    
    % determine offsets (microns) between object bounding box center and the
    % destination image center
    d_XobjOffsetFromImgCenterMicron= d_XObjcenterPtMicron - d_XImg;
    d_YobjOffsetFromImgCenterMicron= d_YObjcenterPtMicron - d_YImg;
    
    
    % determine the pixel coordinates of the center of the bounding box in
    % the DESTINATION image
    
    switch RowIsDim
        case '+y' % should be this for E157
            d_RowOffsetOfObjCenterFromImgCenter=   round(d_YobjOffsetFromImgCenterMicron/dPixelSizeMicrons);
        case '-y'
            d_RowOffsetOfObjCenterFromImgCenter= - round(d_YobjOffsetFromImgCenterMicron/dPixelSizeMicrons);
        otherwise
            error('cannot handle RowIsDim=%s',RowIsDim)
    end
    switch ColIsDim
        case '+x'
            d_ColOffsetOfObjCenterFromImgCenter=   round(d_XobjOffsetFromImgCenterMicron/dPixelSizeMicrons);
        case '-x' % should be this for E157
            d_ColOffsetOfObjCenterFromImgCenter= - round(d_XobjOffsetFromImgCenterMicron/dPixelSizeMicrons);
        otherwise
            error('cannot handle ColIsDim=%s',ColIsDim)
    end
    
    d_objCenterPtRowCol= d_RowColOfImgCenter + [d_RowOffsetOfObjCenterFromImgCenter, d_ColOffsetOfObjCenterFromImgCenter];
    assert(all(mod(d_objCenterPtRowCol,1)==0))
    % calculate the bounding box in pixel coordinates of the DESTINATION image
    dBbox=nan(1,4); % [colStart, rowStart, colWidth, rowHeight]
    dMaskCroppedSize=round(sPixelSizeMicrons/dPixelSizeMicrons*[size(sMaskCropped)]);
    
    dBbox(2)= round(d_objCenterPtRowCol(1) - dMaskCroppedSize(1)/2); %rowStart
    dBbox(1)= round(d_objCenterPtRowCol(2) - dMaskCroppedSize(2)/2); %colStart
    dBbox(4)= dMaskCroppedSize(1); %rowHeight
    dBbox(3)= dMaskCroppedSize(2); %colWidth
    
    %[row, col; row, col; row,col; row,col];
    dCornerPts=nan(4,2);
    %                   ROW                     COL
    dCornerPts(1,:)=[dBbox(2)            ,dBbox(1)            ];     %topleft
    dCornerPts(2,:)=[dBbox(2)+dBbox(4)-1 ,dBbox(1)            ];     %botleft
    dCornerPts(3,:)=[dBbox(2)            ,dBbox(1)+dBbox(3)-1 ];     %topright
    dCornerPts(4,:)=[dBbox(2)+dBbox(4)-1 ,dBbox(1)+dBbox(3)-1 ];     %botright
    
    % determine if the bounding box is within bounds of this image
    rowInRange=all([dCornerPts(:,1)>=1,dCornerPts(:,1)<=dImgSize(1)],2); % rows of corners in range
    colInRange=all([dCornerPts(:,2)>=1,dCornerPts(:,2)<=dImgSize(2)],2); % cols of corners in range
    rowAndColInRange=all([rowInRange,colInRange],2); % each of the 4 corner points
    
    if sum(rowAndColInRange)==4
        %fprintf("the bounding box is INSIDE destination image %i\n",indDestImg)
        % all of the corners of the bounding box are within the destination image
        out_dImageFileSuffix(iObj)=indDestImg;
        out_dBoundingBoxes(iObj,1:4)=dBbox;
        
        % dMask cropped
        dMaskCropped=imresize(sMaskCropped,dMaskCroppedSize);
        out_dMasksCropped(iObj)={dMaskCropped};
        
        % dMask full
        dMask=false(dImgSize);
        dMask(dCornerPts(1,1):dCornerPts(4,1),dCornerPts(1,2):dCornerPts(4,2))=dMaskCropped;
        out_dMasks(iObj)={dMask};
        
        
    else % this segmentation is out of range of even the closest destination image, so fill in variables as nan
        
        
        % for fun
        if sum(rowAndColInRange)>=1
            fprintf("the bounding box is partly outside of the destination image that it's closest to, which is %i\n",indDestImg)
            % bounding box overlaps this image, but it's not fully inside of it
        elseif sum(rowAndColInRange)==0
            display('the bounding box is completely outside of any destination image')
            % bounding box does not overlap with the image at all
        end
    end
    % if the bounding box is out of bounds, then output nan for this
    % object. Otherwise calculate a new cropped mask for this
    
    
    %     clf; figure(1); subplot(1,2,1)
    %     imshow(sMaskCropped)
    %     subplot(1,2,2)
    %     imshow(imresize(sMaskCropped,dMaskCroppedSize))
    %
end


end