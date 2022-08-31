function [pairsForAngle,pairsForFlip,angleDim]=chooseTilePairsForAngleAndRowFlip(XDimConsecPairs,YDimConsecPairs,XDimNonConsecPairs,YDimNonConsecPairs,numTilePairsToTryToFindAngleFor,numTilePairsToTryToFindFlippedDimFor,consecutivePairBiasFactor)
% based on overlapping tile pairs from getEstimatedOverlappingTilePairs,
% decide which ones to use to determine camera angle
%
% preference is for the dimension with the most consecutive pairs, since relative XY positions are
% (maybe) more accurate
%
% however, below the threshold numImagePairsToTryToFindAngleFor consider
% using the other dimension. consecutive pairs are weighted by a bias
% factor consecutivePairBiasFactor, like consecutivePairBiasFactor=1.2
%
%
rng(0)
% we prefer consecutive image pairs, since I trust the relative XY
% positions more however if these are too few then include even
% non-consecutive pairs
if size(XDimConsecPairs,1) > size(YDimConsecPairs,1)
    % then maybe we should use xdim
    if size(XDimConsecPairs,1)<numTilePairsToTryToFindAngleFor
        % but there aren't enough XDimConsecPairs, so now consider YDim.
        % Bias the result towards consec pairs though by a factor of 1.2 (arbitrary)
        if  (consecutivePairBiasFactor*size(XDimConsecPairs,1) + 1*size(XDimNonConsecPairs,1)) > (consecutivePairBiasFactor*size(YDimConsecPairs,1) + 1*size(YDimNonConsecPairs,1))
            % so when weighting for superiority of consecutive pairs, still choose x
            angleDim='x';
            chooseK=min(numTilePairsToTryToFindAngleFor-size(XDimConsecPairs,1),size(XDimNonConsecPairs,1));
            pairsForAngle=[XDimConsecPairs;XDimNonConsecPairs(sort(randperm(size(XDimNonConsecPairs,1),chooseK)),:)];
        else % ydim now looks better, due to greater number of NonConsec pairs
            angleDim='y';
            chooseK=min(numTilePairsToTryToFindAngleFor-size(YDimConsecPairs,1),size(YDimNonConsecPairs,1));
            pairsForAngle=[YDimConsecPairs;YDimNonConsecPairs(sort(randperm(size(YDimNonConsecPairs,1),chooseK)),:)];
        end
    else
        % there are enough image pairs here in XDimConsecPairs
        angleDim='x';
        pairsForAngle=XDimConsecPairs(sort(randperm(size(XDimConsecPairs,1),numTilePairsToTryToFindAngleFor)),:);
    end
else
    % then maybe we should use ydim
    if size(YDimConsecPairs,1)<numTilePairsToTryToFindAngleFor
        % but there aren't enough YDimConsecPairs, so now consider XDim.
        % Bias the result towards consec pairs though by a factor of 1.2 (arbitrary)
        if (consecutivePairBiasFactor*size(YDimConsecPairs,1) + 1*size(YDimNonConsecPairs,1)) > (consecutivePairBiasFactor*size(XDimConsecPairs,1) + 1*size(XDimNonConsecPairs,1))
            % so when weighting for superiority of consecutive pairs, still choose y
            angleDim='y';
            chooseK=min(numTilePairsToTryToFindAngleFor-size(YDimConsecPairs,1),size(YDimNonConsecPairs,1));
            pairsForAngle=[YDimConsecPairs;YDimNonConsecPairs(sort(randperm(size(YDimNonConsecPairs,1),chooseK)),:)];
        else
            % xdim now looks better, due to greater number of NonConsec pairs
            angleDim='x';
            chooseK=min(numTilePairsToTryToFindAngleFor-size(XDimConsecPairs,1),size(XDimNonConsecPairs,1));
            pairsForAngle=[XDimConsecPairs;XDimNonConsecPairs(sort(randperm(size(XDimNonConsecPairs,1),chooseK)),:)];
        end
    else
        % there are enough image pairs here in YDimConsecPairs
        angleDim='y';
        pairsForAngle=YDimConsecPairs(sort(randperm(size(YDimConsecPairs,1),numTilePairsToTryToFindAngleFor)),:);
    end
end

%% now get other pairs for determining whether camera is flipped
if strcmp(angleDim,'x')
    if size(YDimConsecPairs,1)<numTilePairsToTryToFindFlippedDimFor
        % then use some nonconsec Y pairs
        chooseK=min(numTilePairsToTryToFindFlippedDimFor-size(YDimConsecPairs,1),size(YDimNonConsecPairs,1));
        pairsForFlip=[YDimConsecPairs;YDimNonConsecPairs(sort(randperm(size(YDimNonConsecPairs,1),chooseK)),:)];
    else
        % have enough consec Y pairs
        pairsForFlip=YDimConsecPairs(sort(randperm(size(YDimConsecPairs,1),numTilePairsToTryToFindFlippedDimFor)),:);
    end
elseif strcmp(angleDim,'y')
    if size(XDimConsecPairs,1)<numTilePairsToTryToFindFlippedDimFor
        % then use all Xdim pairs
        chooseK=min(numTilePairsToTryToFindFlippedDimFor-size(YDimConsecPairs,1),size(XDimNonConsecPairs,1));
        pairsForFlip=[XDimConsecPairs;XDimNonConsecPairs(sort(randperm(size(XDimNonConsecPairs,1),chooseK)),:)];
    else
        % have enough consec X pairs
        pairsForFlip=XDimConsecPairs(sort(randperm(size(XDimConsecPairs,1),numTilePairsToTryToFindFlippedDimFor)),:);
    end
else
    error('error in internal logic')
end


end