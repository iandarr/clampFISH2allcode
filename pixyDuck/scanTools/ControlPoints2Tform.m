function tform=ControlPoints2Tform(ControlPoints)
% 
numCpPairs=length(ControlPoints);
%fixed_rstage_centerXY=reshape([ControlPoints(:).fixed_rstage_centerXY],2,numCpPairs)';
move_stage_centerXY=reshape([ControlPoints(:).move_stage_centerXY],2,numCpPairs)';
move_OffsetXY=double(reshape([ControlPoints(:).move_OffsetXY],2,numCpPairs)');

%tform=fitgeotrans(move_stage_centerXY,[move_stage_centerXY+move_OffsetXY],'nonreflectivesimilarity'); % this ordering, the inverse trasnform can be used: registeredXY=transformPointsInverse(tform,unregisteredXY)
%tform=fitgeotrans(move_stage_centerXY,[move_stage_centerXY+move_OffsetXY],'pwl'); % pwl= PiecewiseLinearTransformation2D with this, we only have transformPointsInverse which can go to fixed-->move so for me ordering is (fixed,moving) which is reverse of how data is handled internally
[tform,inliersVect]=estimateGeometricTransform2D(move_stage_centerXY,[move_stage_centerXY+move_OffsetXY],'rigid','MaxDistance',25); % note: 'MaxDistance' is in um since that's stage coordinate units. For reference my nuclei are about 10-20um in diameter
if sum(inliersVect)~=numCpPairs % all should be inliers. 
   warning('only have %i inliers (out of %i control poitns) in registration', sum(inliersVect),numCpPairs)
end
% could also try these transforms:
% PolynomialTransformation2D
% PiecewiseLinearTransformation2D

end