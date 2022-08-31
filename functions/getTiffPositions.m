function [Xpositions,Ypositions]=getTiffPositions(tiffFileList)
% gets positions of Tiff files from Tiff exported from NIS-Elements
% 
% Xpositions and Ypositions are column vectors

if ischar(tiffFileList)
    tiffFileList=cellstr(tiffFileList);
end


nTiffFiles=length(tiffFileList);
Xpositions=nan(nTiffFiles,1);
Ypositions=nan(nTiffFiles,1);


for i=1:nTiffFiles
	
    filename=tiffFileList{i};
    
    if ~isfile(filename)
       error("could not find file %s",filename) 
    end
    %obj = Tiff(filename);
    %obj.getTagNames()
    %getTag(obj,'XPosition')
    %getTag(obj,obj.TagID.XPosition)
    
    info=imfinfo(filename);
    ImageDescription=info(1).ImageDescription;
   
   % X position
   positionXexpression='PositionX="';
   positionXcell=regexp(ImageDescription,[positionXexpression,'[0-9,.,-]*'],'match');
   if length(positionXcell)~=1
        error("regexp found more than one PositionX="" in the ImageDescription")
   end
   positionX=str2double(replace(positionXcell{1},positionXexpression,''));
   Xpositions(i)=positionX;
   
   % Y position
   positionYexpression='PositionY="';
   positionYcell=regexp(ImageDescription,[positionYexpression,'[0-9,.,-]*'],'match');
   if length(positionYcell)~=1
        error("regexp found more than one PositionY="" in the ImageDescription")
   end
   positionY=str2double(replace(positionYcell{1},positionYexpression,''));
   Ypositions(i)=positionY;
   
end
end