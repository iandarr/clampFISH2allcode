function middleAngle=middleAngle(AnglesInDegrees)
% middle angle finds the 'middle' angle
% like a median, however if you just take the median you can get wrap-around effects. Like
% for example if you had this list of AnglesInDegrees
%    AnglesInDegrees=[-182 -181 0 181 182]'
% and just took median(AnglesInDegrees) then you'd get zero degrees, while
% what you really want is 180 or -180 degrees
% this is not always solved by wrapping the angles to, say, positive numbers
% because the same thing can happen. For example:
%   AnglesInDegrees=[1 2 180 358 359] where median(AnglesInDegrees) = 180,
%   while what you really want is 0 degrees
%
%
% If there are an even number of angles it will likely return the 
AnglesInDegrees=wrapTo360(AnglesInDegrees(~isnan(AnglesInDegrees))); % can't be nan, wrap to 360
assert(isvector(AnglesInDegrees))
if size(AnglesInDegrees,1)==1 % then input is a row vector (unless scalar)
    AnglesInDegrees=AnglesInDegrees'; % make it column vector
end
numAngles=length(AnglesInDegrees);

AnglesInRadians=wrapTo2Pi(AnglesInDegrees*pi()/180);

% example:
% angleDifference=angdiff(80*pi()/180,90*pi()/180)*180/pi()
% angleDifference = 10
angDiffsMat=angdiff(repmat(AnglesInRadians,1,length(AnglesInRadians)),repmat(AnglesInRadians',length(AnglesInRadians),1));

%keepDiffMat=ones(numAngles);
%keepDiffMat(logical(tril(ones(numAngles))))=nan; % nan on the diagonal (self) and in lower-left triangle (duplicates)
keepDiffMat=ones(numAngles);
keepDiffMat(logical(diag(ones(1,numAngles))))=nan; % nan along diagonal
angDiffsMat2=keepDiffMat.*angDiffsMat; % now self is nan
angDiffsAbsSum=nansum(abs(angDiffsMat2),2); % sum of absolute values of angular differences;
[minVal,~]=min(angDiffsAbsSum);
idxMin=angDiffsAbsSum==minVal; % do it this way in case there are 2+ angles with min value

idxMin=abs(angDiffsAbsSum-minVal) < 1e4*eps(min(min(abs(angDiffsAbsSum),abs(minVal)))); % within 1e4 * eps


if sum(idxMin)>1
    %fprintf('>1 min angle...\n')
    if sum(idxMin)==2
        %fprintf(' 2 angles to find midpoint of from\n')
        % choose in the middle of the two
        middleAnglesToChooseFromInRadians=AnglesInRadians(idxMin);
        angDiffChosenAngles=angdiff(middleAnglesToChooseFromInRadians(1),middleAnglesToChooseFromInRadians(2));
        middleAngle=180/pi()*(middleAnglesToChooseFromInRadians(1) + angDiffChosenAngles/2);
    else  % more than 2 angles 
        %fprintf(' >2 angles to find midpoint of from\n')
        if all(abs(AnglesInRadians(idxMin)-mean(AnglesInRadians(idxMin)))<=1e4*eps(min(min(abs(AnglesInRadians),abs(mean(AnglesInRadians(idxMin))))))) % could have identical angles?
            middleAngle=180/pi()*AnglesInRadians(find(idxMin,1,'first')); % just return the first of these
            %fprintf(' they are all the same though, return the first\n')
        else
            % cannot handle this case. Could be radially symmetric
            warning("more than 2 angles were found to be closest to all the others. Don't know how to handle this. Will return nan. Can happen when angles are radially symmetric")
            middleAngle=nan;
        end
    end
else
    middleAngle=180/pi()*AnglesInRadians(idxMin);
end


middleAngle=wrapTo180(middleAngle);

end