function fh=packSubplotsTogether(fh,xfactor,yfactor)
    % factor should be greater than 1 to make bigger in that direction
    
for iChild=1:length(fh.Children)
    
    Position=fh.Children(iChild).Position; %initial position
    
    width=Position(3);
    height=Position(4);
    
    centerX=Position(1)+ 0.5 * width;
    centerY=Position(2)+ 0.5 * height;
    
    Position(1)=centerX - 0.5 * width * xfactor;
    Position(2)=centerY - 0.5 * height *yfactor;
    
    Position(3)= width*xfactor;
    Position(4)= height*xfactor;
    
    fh.Children(iChild).Position=Position;
end

end