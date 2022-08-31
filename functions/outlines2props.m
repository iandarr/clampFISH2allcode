function Tprops=outlines2props(cellposeOutlinesFilepath,UmPerPx)
% Tprops has columns:
%
assert(isscalar(UmPerPx))

%% read cellpose's outlines .txt file to get polyshape vector. This takes a long time.
polygonVector=d2utils.parseCellposeOutlines(outlinesPath, 'flip',false);

%% find shape properties
% bounding box
[collim,rowlim]=arrayfun(@(x) boundingbox(x),polygonVector,'UniformOutput',false);
collim=cell2mat(collim'); rowlim=cell2mat(rowlim'); % to numeric arrays
bbox=[rowlim(:,1), collim(:,1), rowlim(:,2)-rowlim(:,1)+1, collim(:,2)-collim(:,1)+1]; %[rowstart colstart height width]

% centroid
%[colC,rowC]=arrayfun(@(x) centroid(x),polygonVector,'UniformOutput',false);
[centCol,centRow]=arrayfun(@(x) centroid(x),polygonVector);
centCol=centCol'; centRow=centRow';

% area
areaPoly=arrayfun(@(x) area(x),polygonVector)*(UmPerPx^2);

% perimeter
perim=arrayfun(@(x) perimeter(x),polygonVector)*UmPerPx;

Tprops=table('Size',[length(polygonVector) 6]);
Tprops.centRow=centRow;
Tprops.centCol=centCol;
Tprops.bbox1=bbox(:,1);
Tprops.bbox2=bbox(:,2);
Tprops.bbox3=bbox(:,3);
Tprops.bbox4=bbox(:,4);
Tprops.area_Um2=areaPoly;
Tprops.perimeter_Um=perim;


end