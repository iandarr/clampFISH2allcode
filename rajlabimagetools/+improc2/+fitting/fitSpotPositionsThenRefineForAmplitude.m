function [fittedSpots, fittedBackgLevels] = fitSpotPositionsThenRefineForAmplitude(...
        imgForPositionDetermination, imgForAmplitudeFit, ...
        spotXGuesses, spotYGuesses, spotZPlanes)
    
    
    numSpots = length(spotXGuesses);
    
    stageOneFitter = improc2.fitting.GaussianSpotFitter(imgForPositionDetermination);
    stageOneFitter.halfLengthOfRegionToFit = 4;
    stageOneFitter.xMaxDeviationFromGuess = 2;
    stageOneFitter.yMaxDeviationFromGuess = 2;
    stageOneFitter.sigmaGuess = 1.5;
    stageOneFitter.sigmaLimits = [0 3.5];
    stageOneFitter.offsetLimits = [0 Inf];
    
    fittedSpotsStageOne = preAllocateGaussian2dSpotArray(numSpots);
    
    for i = 1:numSpots
        stageOneFitter.xGuess = spotXGuesses(i);
        stageOneFitter.yGuess = spotYGuesses(i);
        stageOneFitter.spotZPlane = spotZPlanes(i);
        fittedSpotsStageOne(i) = stageOneFitter.fitSpot();
    end

    stageTwoFitter = improc2.fitting.GaussianSpotFitter(imgForAmplitudeFit);
    stageTwoFitter.halfLengthOfRegionToFit = 7;
    stageTwoFitter.xMaxDeviationFromGuess = 0.2;
    stageTwoFitter.yMaxDeviationFromGuess = 0.2;
    stageTwoFitter.sigmaGuess = 1.5;
    stageTwoFitter.sigmaLimits = [0 3.5];
    stageTwoFitter.offsetLimits = [0 Inf];
    
    fittedSpots = preAllocateGaussian2dSpotArray(numSpots);
    fittedBackgLevels = zeros(1, numSpots);
    
    for i = 1:numSpots
        stageTwoFitter.xGuess = fittedSpotsStageOne(i).xCenter;
        stageTwoFitter.yGuess = fittedSpotsStageOne(i).yCenter;
        stageTwoFitter.spotZPlane = spotZPlanes(i);
        [fittedSpots(i), fittedBackgLevels(i)] = stageTwoFitter.fitSpot();
    end
end

function spotArray = preAllocateGaussian2dSpotArray(numSpots)
    spotArray(1, numSpots) = improc2.fitting.Gaussian2dSpot();
end