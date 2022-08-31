function flag=shiftCurrentObjSegmentationMask(objectHandle,xShift,yShift)
% shiftCurrentObjSegmentationMask(xShift,yShift)
%
% takes the segmentation mask of the current object and shifts it by iShift
% pixels down and jShift pixels right
%
% objectHandle should be gotten like this:
%   tools=improc2.launchImageObjectTools;
%   tools.iterator.goToFirstObject() <-- or whatever object you want
%   objectHandle=tools.objectHandle;

ImageObjectBaseData=objectHandle.getData('imageObject','improc2.dataNodes.ImageObjectBaseData');
imfilemaskOrig=ImageObjectBaseData.imfilemask;


[maxIndY,maxIndX]=size(imfilemaskOrig);
%imshow(imfilemaskOrig)

[indY,indX]=find(imfilemaskOrig);

newIndX=indX+xShift;
newIndY=indY+yShift;

imfilemaskNew=false(size(imfilemaskOrig));

segmentationIsInBounds=all(all([newIndX>=1,newIndY>=1,newIndX<=maxIndX,newIndY<=maxIndY],2),1);

if segmentationIsInBounds
    newIndLinear=sub2ind(size(imfilemaskOrig),newIndY,newIndX);
    imfilemaskNew(newIndLinear)=true;
    flag=1; %success
else % not in bounds - can't use this 
    
    imfilemaskNew(1:10,1:10)=true; % just 1x1 pixel segmentation
    flag=0; %fail
end
    
%imshow(imfilemaskNew)

%imshowpair(imfilemaskOrig,imfilemaskNew,'montage')

%imshowpair(imfilemaskOrig, imfilemaskNew,'Scaling','joint');

ImageObjectBaseDataNew=ImageObjectBaseData;
ImageObjectBaseDataNew.imfilemask=imfilemaskNew;
objectHandle.setData(ImageObjectBaseDataNew,'imageObject','improc2.dataNodes.ImageObjectBaseData')



end