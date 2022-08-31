function [fittedSpots, fittedBackgLevels] = fitSpotPositionsThenRefineForAmplitude2(...
        imgForPositionDetermination, imgForAmplitudeFit, ...
        spotXGuesses, spotYGuesses, spotZPlanes)
    % copied from rajlabimagetools' fitSpotPositionsThenRefineForAmplitude
    % with parameters modified
    
    numSpots = length(spotXGuesses);
    
    stageOneFitter = improc2.fitting.GaussianSpotFitter(imgForPositionDetermination);
    % GaussianSpotFitter2 may have edited parameters from dentist2 code
    stageOneFitter.halfLengthOfRegionToFit = 2;
    stageOneFitter.xMaxDeviationFromGuess = 1;
    stageOneFitter.yMaxDeviationFromGuess = 1;
    stageOneFitter.sigmaGuess = 0.65;
    stageOneFitter.sigmaLimits = [0.65 0.75];
    stageOneFitter.offsetLimits = [0 Inf];
    
    fittedSpotsStageOne = preAllocateGaussian2dSpotArray(numSpots);
    
    for i = 1:numSpots
        stageOneFitter.xGuess = spotXGuesses(i);
        stageOneFitter.yGuess = spotYGuesses(i);
        stageOneFitter.spotZPlane = spotZPlanes(i);
        fittedSpotsStageOne(i) = stageOneFitter.fitSpot();
    end

    stageTwoFitter = improc2.fitting.GaussianSpotFitter(imgForAmplitudeFit);
    % GaussianSpotFitter2 may have edited parameters from dentist2 code
    stageTwoFitter.halfLengthOfRegionToFit = 2;
    stageTwoFitter.xMaxDeviationFromGuess = 0.2;
    stageTwoFitter.yMaxDeviationFromGuess = 0.2;
    stageTwoFitter.sigmaGuess = 0.65;
    stageTwoFitter.sigmaLimits = [0.65 0.75];
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