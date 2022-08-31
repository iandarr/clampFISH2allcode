function SURFPointsInCenter=returnOnlyCentralSURFPoints(SURFPoints,imgHeight,imgWidth,linearFractionOfImage)


leftmostPx=round(imgWidth*(1-linearFractionOfImage)/2);
rightmostPx=round(imgWidth*(1-(1-linearFractionOfImage)/2));

topmostPx=round(imgHeight*(1-linearFractionOfImage)/2);
bottommostPx=round(imgHeight*(1-(1-linearFractionOfImage)/2));

XYpos=[SURFPoints(:).Location];
X=XYpos(:,1);
Y=XYpos(:,2);

inCentralRegion=all([leftmostPx<=X,rightmostPx>=X,topmostPx<=Y,bottommostPx>=Y],2);

SURFPointsInCenter=SURFPoints(inCentralRegion);
end