function polyvec=polyshapeVector(X_nx4,Y_nx4)
% outputs a Nx1 array of polyshapes
%
% each row of X_nxy and Y_nx4 are 1x4 arrays that define the x and y
% locations of the corners of the polyshapes

numPoly=size(X_nx4,1);
assert(size(X_nx4,2)==4)
assert(isequal(size(X_nx4),size(Y_nx4)))

polyvec=polyshape.empty([0,numPoly]);
for i=1:numPoly
    polyvec(i)=polyshape(X_nx4(i,1:4),Y_nx4(i,1:4));
end

end