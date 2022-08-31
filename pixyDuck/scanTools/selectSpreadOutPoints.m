function indSelected=selectSpreadOutPoints(XYall,numAdditional,indSelected)
% given XY locations of points (XYall) and the indices of existing selected
% points (indSelected), select an additional numAdditional points iteratively
% such that each successive point is the farthest away from already
% selected points
% 
% indSelected should be a column vector
% 
assert(size(XYall,2)==2)
assert(size(indSelected,2)==1)


numSelected=size(indSelected,1);
assert( (numAdditional+numSelected) < size(XYall,1)) % there must be enough points to choose from

% recursive
if numAdditional>0
    
    XYselected=XYall(indSelected,:);
    [distMaxMin,indNew]=max(min(pdist2(XYall,XYselected),[],2));
    
    
    indSelected=[indSelected;indNew];
    numAdditional=numAdditional-1;
    
    indSelected=selectSpreadOutPoints(XYall,numAdditional,indSelected);
end


end