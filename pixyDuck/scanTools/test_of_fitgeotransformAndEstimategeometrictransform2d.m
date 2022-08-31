


clc

XYorig=[0 1; 5 8; 10 1];
XYnew=[0 1; 5 8; 10 1].*[1 -1]; % flip in Y

fprintf('--- estimateGeometricTransform2D ridig (no reflection)---\n')
[tform2dFn1,inlierIndex]=estimateGeometricTransform2D(XYorig,XYnew,'rigid')  % similarity does not allow reflection
tform2dFn_det=det(tform2dFn1.T)
fprintf('\n\n')

fprintf('--- estimateGeometricTransform2D similarity (no reflection)---\n')
[tform2dFn2,inlierIndex]=estimateGeometricTransform2D(XYorig,XYnew,'similarity')  % similarity does not allow reflection
tform2dFn_det=det(tform2dFn2.T)
fprintf('\n\n')

fprintf('--- estimateGeometricTransform2D affine (allows reflection)---\n')
[tform2dFn3,inlierIndex]=estimateGeometricTransform2D(XYorig,XYnew,'affine')  % affine does allow reflection
tform2dFn_det=det(tform2dFn3.T)
fprintf('\n\n')

fprintf('--- fitgeotrans similarity (allows reflection) ---\n')
[tformGeo]=fitgeotrans(XYorig,XYnew,'similarity') % similarity allows reflection
tformGeo_det=det(tformGeo.T)
fprintf('\n\n')

fprintf('--- fitgeotrans nonreflectivesimilarity (doesnt allow reflection) ---\n')
[tformGeoNonReflective]=fitgeotrans(XYorig,XYnew,'nonreflectivesimilarity') % now we suppress reflection
tformGeoNonReflective_det=det(tformGeoNonReflective.T)


