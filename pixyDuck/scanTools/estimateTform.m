function [tform,numInliers,statusFlag]=estimateTform(img_move,img_fixed,varargin)
% tform=estimateTform(img_move,img_fixed)
%
% tform=estimateTform(img_move,img_fixed,'Name',Value,...)
%
% this is a wrapper function for MATLAB's detectSURFFeatures, extractFeatures, matchFeatures, and estimateGeometricTransform2D
%
% 
% Name-Value Pairs
% 
%  'transformType'
%   'rigid' (default) |'similarity'|'affine'|'projective' | cell array with these values for multiple types of tform
%       Example: transformType='rigid'      (translation and rotation; no reflection)
%       Example: transformType='similarity' (translation, rotation, and scaling; no reflection)
%       Example: transformType='affine'     (translation, rotation, scaling, shearing; reflection is allowed)
%       Example: transformType={'rigid'; 'similarity';'affine'}
%           in this case, output will be:
%               tform        3x1 cell array
%               numInliers   3x1 numeric array
%               statusFlag   3x1 numeric array
%
%  'ROI_img_fixed' - detectSURFFeatures will only find features in this fixed image ROI
%    1x4 numeric vector in format [col1 row1 width height]
%
%  'ROI_img_fixed' - detectSURFFeatures will only find features in this move image ROI
%    1x4 numeric vector in format [col1 row1 width height]
%   
%   'minNumberOfInliers' - if numInliers is below this threshold, this function will return NaN values with statusFlag = 3
%    numeric scalar (default 2)
%
%   'DisplayMatchingPointPairsUsedInTform' - will plot the features found and how they have been matched by estimateGeometricTransform2D
%    true | false
%
%   ----------- Inputs related to image scaling ----------
%   'intensityNormType' - pixel intensities are multiplied by a factor such that the N'th percentile is equal to one. Above 1 is set equal to 1.
%   'none' (default) | 'SameScaling' | 'EachImage'
%   
%   'intensityNormHighPercentile'
%    scalar numeric (default 99.999)
%
%   Example:
%       estimateTform(img_move,img_fixed,'IntensityNormType','SameScaling','intensityNormHighPercentile',99.999)
%       will take
%       normFactor=1/max(prctile(img_fixed,99.999,'all'),prctile(img_move,99.999,'all'));
%   
%       and then scale them like this:
%           imgFixed=normFactor*imgFixed;
%           imgMove=normFactor*imgMove;
%   
%       whereas 'SameScaling' would use each individual image's 99.999th
%       percentile in their own scalefactor
%
%  -------- inputs that get passed to matchFeatures --------
%  'MatchThreshold'   A scalar T, 0 < T <= 100, that specifies the
%                      distance threshold required for a match. A pair of
%                      features are not matched if the distance between
%                      them is more than T percent from a perfect match.
%                      Increase T to return more matches.
%
%                      matchFeatures default (for non-binary features) is 1
%
%
%   'MaxRatio'         A scalar R, 0 < R <= 1, specifying a ratio threshold
%                      for rejecting ambiguous matches. Increase R to
%                      return more matches.
%
%                      Default: 0.6
% 
% ------- statusFlag meanings ---------
%    statusFlag set from estimateGeometricTransform2D:
%           0: No error.
%           1: matchedPoints1 and matchedPoints2 do not contain enough points.
%           2: Not enough inliers have been found (per estimateGeometricTransform2D)
%
% Example 2:
%   MetricThreshold=30;
%   topStrongestSURFPoints=100;
%   linearFractionOfImage=0.6;
%   DisplayFeaturePoints=true;
%   DisplayMatchingPointPairsUsedInTform=true;
%   tform=estimateTform(img_move,img_ref,'MetricThreshold',MetricThreshold,'UseCentralFractionOfImage',linearFractionOfImage,'NumberOfStrongestSURFFeaturesToUse',topStrongestSURFPoints);
%

ip=inputParser;
ip.KeepUnmatched=true;
%ip.addParameter('MetricThreshold',[],@isnumeric)
ip.addParameter('NumberOfStrongestSURFFeaturesToUse',[],@isnumeric)
%ip.addParameter('DisplayFeaturePoints',false,@islogical)
ip.addParameter('DisplayMatchingPointPairsUsedInTform',false,@islogical)
ip.addParameter('transformType','rigid',@(x) all(ismember(x,{'rigid','similarity','affine','projective'})))
ip.addParameter('featuresMatchThreshold',5,@(x) all([isnumeric(x),x>0,x<=100])) % passed to matchFeatures, whose numeric default is 1. 
ip.addParameter('featuresMatchMaxRatio',0.6,@(x) and(x>0,x<=1))
ip.addParameter('MaxNumTrials',2000,@(x) x>=500)
ip.addParameter('ROI_img_fixed',[],@(x) and(isnumeric(x),isequal(size(x),[1 4])))
ip.addParameter('ROI_img_move', [],@(x) and(isnumeric(x),isequal(size(x),[1 4])))

ip.addParameter('intensityNormType','none',@(x) any(strcmp(x,{'none','SameScaling','EachImage'})))
ip.addParameter('intensityNormHighPercentile',99.999,@(x) and(isnumeric(x),isscalar(x)))
%ip.addParameter('ROI_img_fixed',[])
%ip.addParameter('ROI_img_move', [])

% handle empty varargin case
for i=1:numel(varargin)
    if  isempty(varargin{i})
        keepIdx=true(1,numel(varargin)); keepIdx(i)=false; varargin=varargin(keepIdx);
    end
end

% varginPassthrough is for detectSURFFeatures
ip.parse(varargin{:});
% handle empty varargin case
if isempty(varargin)
    % skip input parsing
    varginPassthrough=struct();
else %something was in varargin, pass on what wasn't used
    varginPassthrough=ip.Unmatched;    
end
if ~ismember('MetricThreshold',fields(varginPassthrough))
    varginPassthrough.MetricThreshold=50; % default
end

%MetricThreshold=ip.Results.MetricThreshold;
NumberOfStrongestSURFFeaturesToUse=ip.Results.NumberOfStrongestSURFFeaturesToUse;
%DisplayFeaturePoints=ip.Results.DisplayFeaturePoints;
DisplayMatchingPointPairsUsedInTform=ip.Results.DisplayMatchingPointPairsUsedInTform;
transformType=ip.Results.transformType;
ROI_img_fixed=ip.Results.ROI_img_fixed;
ROI_img_move=ip.Results.ROI_img_move;
featuresMatchThreshold=ip.Results.featuresMatchThreshold;
featuresMatchMaxRatio=ip.Results.featuresMatchMaxRatio;
intensityNormType=ip.Results.intensityNormType;
intensityNormHighPercentile=ip.Results.intensityNormHighPercentile;
MaxNumTrials=ip.Results.MaxNumTrials; %for estimateGeometricTransform2D

%% ROIs are specified
% fixed
varginPassthrough_img_fixed=varginPassthrough;
if ~isempty(ROI_img_fixed)
    varginPassthrough_img_fixed.ROI=ROI_img_fixed;
end
% move
varginPassthrough_img_move=varginPassthrough;
if ~isempty(ROI_img_move)
   varginPassthrough_img_move.ROI=ROI_img_move;
end

%% scale image intensity
if strcmp(intensityNormType,'none')
    % do nothing
else
    intensityHigh_fixed=prctile(img_fixed,intensityNormHighPercentile,'all');
    intensityHigh_move=prctile(img_move,intensityNormHighPercentile,'all');
    
    if strcmp(intensityNormType,'SameScaling')
        scaleFactorMove=1/max(intensityHigh_fixed,intensityHigh_move);
        scaleFactorFixed=scaleFactorMove;
    elseif strcmp(intensityNormType,'EachImage')
        scaleFactorMove=1/intensityHigh_move;
        scaleFactorFixed=1/intensityHigh_fixed;
    else
        error('error in logic')
    end
    img_move=scaleFactorMove*img_move;
    img_move(img_move>1)=1;
    img_fixed=scaleFactorFixed*img_fixed;
    img_fixed(img_fixed>1)=1;   
end

%% get features with detectSURFFeatures
% this is the most time-consuming step

if isempty(varginPassthrough_img_fixed)
    featurePts_fixed= detectSURFFeatures(img_fixed); % default is MetricThreshold=1000
else
    featurePts_fixed= detectSURFFeatures(img_fixed,varginPassthrough_img_fixed);
end

if isempty(varginPassthrough_img_move)
    featurePts_move= detectSURFFeatures(img_move);% default is MetricThreshold=1000
else
    featurePts_move= detectSURFFeatures(img_move,varginPassthrough_img_move);% default is MetricThreshold=1000
end

num_featurePts_fixed_found=length(featurePts_fixed);
num_featurePts_move_found=length(featurePts_move);

%% Select only strongest features for tform
if isempty(NumberOfStrongestSURFFeaturesToUse)
    % do nothing - use all the features
else
    featurePts_fixed=featurePts_fixed.selectStrongest(NumberOfStrongestSURFFeaturesToUse);
    featurePts_move=featurePts_move.selectStrongest(NumberOfStrongestSURFFeaturesToUse);
end

% %% optionally show the features
% if DisplayFeaturePoints
%     
%     figure('Position', [10 10 1000 500])
%     subplot(1,2,1);
%     hold on;
%     set(gca, 'YDir','reverse')
%     set(gca, 'Position',[0.05 0.1 0.4 0.95])
%     axis equal
%     title('ref image, features')
%     imagesc(img_fixed)
%     
%     plot(featurePts_ref)
%     
%     subplot(1,2,2);
%     hold on;
%     set(gca, 'YDir','reverse')
%     set(gca, 'Position',[0.55 0.1 0.4 0.95])
%     axis equal
%     title('move image, features')
%     imagesc(img_move)
%     plot(featurePts_move)
% end



%% get point descriptors

[features_ref,  validPts_ref]  = extractFeatures(img_fixed,  featurePts_fixed);
[features_move, validPts_move] = extractFeatures(img_move, featurePts_move);

%% match the features

indexPairs = matchFeatures(features_ref, features_move,'MatchThreshold',featuresMatchThreshold,'MaxRatio',featuresMatchMaxRatio);

numMatchedFeatures=size(indexPairs,1);

% retrieve locations of corresponding points for each image
matched_ref  = validPts_ref(indexPairs(:,1));
matched_move = validPts_move(indexPairs(:,2));


%% optionally show putative matched features (but this is before outlier removal)
% showMatchedFeatures(img_ref,img_move,matched_ref,matched_move);


%% estimate transformation using the statistically robust M-estimator SAmple Consensus (MSAC) algorithm

if ischar(transformType)
    transformType={transformType};
end
numTransformTypes=length(transformType);

tform=cell(numTransformTypes,1);
inlierIdx=cell(numTransformTypes,1);
statusFlag=nan(numTransformTypes,1);
numInliers=nan(numTransformTypes,1);

for iTransformType=1:numTransformTypes
    
    transformType_this=transformType{iTransformType};
    
    % call estimateGEometricTransform2D
    warning('off','vision:ransac:maxTrialsReached')
    [tform_this, inlierIdx_this,statusFlag_this] = estimateGeometricTransform2D(...
        matched_move, matched_ref, transformType_this,'MaxNumTrials',MaxNumTrials);
    warning('on','vision:ransac:maxTrialsReached')
    
    tform{iTransformType}=tform_this;
    inlierIdx{iTransformType}=inlierIdx_this;
    numInliers(iTransformType)=sum(inlierIdx_this,'all');
    statusFlag(iTransformType)=statusFlag_this;
end

if numTransformTypes==1
    tform=tform{1}; % just output the transform object, not in a cell
    inlierIdx=inlierIdx{1}; % just output the transform object, not in a cell
end
%% optionally display inliers used in the computation of the transformation

if DisplayMatchingPointPairsUsedInTform
    

    axesSide=0.3;
    lowCLimPercentile=5;
    highCLimPercentile=99.99;
    % ax1 fixed image, features
    fh=figure;
    figureHeight=500;
    set(fh,'Position',[500 500 figureHeight*3 figureHeight])

    %ax1=subplot(1,3,1); %hold on; axis equal;
    %set(ax1,'Position',[0.1 0.1 axesSide axesSide*3],             'YDir','reverse','PositionConstraint','innerposition');hold on; axis equal;
    %ax2=subplot(1,3,2); %hold on; axis equal;
    %set(ax2,'Position',[0.15+axesSide 0.1 axesSide axesSide*3],    'YDir','reverse','PositionConstraint','innerposition');hold on; axis equal;
    %ax3=subplot(1,3,3); %hold on; axis equal;
    %set(ax3,'Position',[0.20+2*axesSide 0.1 axesSide axesSide*3], 'YDir','reverse','PositionConstraint','innerposition');hold on; axis equal;
    
    %ax1=axes(fh,'Position',[0.05 0.05 axesSide axesSide*3],'YDir','reverse','PositionConstraint','innerposition'); hold on; axis equal;
    %ax2=axes(fh,'Position',[0.10+axesSide 0.05 axesSide axesSide*3],'YDir','reverse','PositionConstraint','innerposition'); hold on; axis equal;
    %ax3=axes(fh,'Position',[0.15+2*axesSide 0.05 axesSide axesSide*3],'YDir','reverse','PositionConstraint','innerposition'); hold on; axis equal;
    %tlh=tiledlayout(1,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
    tlh=tiledlayout(1,2+numTransformTypes, 'Padding', 'compact', 'TileSpacing', 'compact'); 
    
    % img fixed, with all features
    ax(1)=nexttile(tlh);
    axis equal; set(gca,'YDir','reverse'); hold on;
    %axes(ax1);
    cmapRed=[1 0 0].*linspace(0,1,256)';
    titleStr={'fixed image, features';...
            sprintf('MetricThreshold=%g',varginPassthrough.MetricThreshold);...
             sprintf('features: %i found, %i used, %i matched',num_featurePts_fixed_found,length(featurePts_fixed),numMatchedFeatures);...
             sprintf('intensityNormType=%s',intensityNormType)};
    title(titleStr)
    DisplayRange=[prctile(img_fixed,lowCLimPercentile,'all'), prctile(img_fixed,highCLimPercentile,'all')];
    h1=imshow(img_fixed,'Parent',ax(1),'DisplayRange',DisplayRange,'Colormap',cmapRed); %imagesc(ax1,img_ref)
    ax1_XLim=ax(1).XLim; ax1_YLim=ax(1).YLim;
    axes(ax(1)); pts1=plot(featurePts_fixed);
    set(findobj(pts1.Children,'Type','Line'),'Color',[1 0 0]) %red
    ax(1).XLim=ax1_XLim; ax(1).YLim=ax1_YLim;
    ax(1).Visible='on'; ax(1).Color='none';

    
    
    % img move, with all features
    ax(2)=nexttile(tlh);
    %axes(ax2);
    axis equal; set(gca,'YDir','reverse'); hold on;
    cmapCyan=[0 1 1].*linspace(0,1,256)';
    titleStr={'move image, features';...
             sprintf('MetricThreshold=%g',varginPassthrough.MetricThreshold);...
             sprintf('features: %i found, %i used, %i matched',num_featurePts_move_found,length(featurePts_move),numMatchedFeatures);...
             sprintf('intensityNormType=%s',intensityNormType)};
    title(titleStr)
    DisplayRange=[prctile(img_move,lowCLimPercentile,'all'), prctile(img_move,highCLimPercentile,'all')];
    h2=imshow(img_move,'Parent',ax(2),'DisplayRange',DisplayRange,'Colormap',cmapCyan); %imagesc(ax2,img_move)
    ax2_XLim=ax(2).XLim; ax2_YLim=ax(2).YLim;
    axes(ax(2)); pts2=plot(featurePts_move);
    set(findobj(pts2.Children,'Type','Line'),'Color',[0 1 0]) %green
    ax(2).XLim=ax2_XLim; ax(2).YLim=ax2_YLim;
    ax(2).Visible='on'; ax(2).Color='none';
    
    
    for iTransformType=1:numTransformTypes
        if iscell(inlierIdx)
            inlierIdx_thisTransformType=inlierIdx{iTransformType};
        else
            inlierIdx_thisTransformType=inlierIdx;
        end
        inlier_move= matched_move(inlierIdx_thisTransformType, :);
        inlier_ref= matched_ref(inlierIdx_thisTransformType, :);
        
        
        % img overlay, with matching features
        axNum=2+iTransformType;
        ax(axNum)=nexttile(tlh);
        axis equal; set(gca,'YDir','reverse'); hold on;
        showMatchedFeatures(img_fixed,img_move,inlier_ref,inlier_move,'Parent',ax(axNum));
        titleStr={  sprintf('%i matched inliers from transformType %s',numInliers(iTransformType),char(transformType(iTransformType)));...
                    sprintf('featuresMatchThreshold=%0.2g',featuresMatchThreshold);...
                    sprintf('featuresMatchMaxRatio=%0.2g',featuresMatchMaxRatio)};
        title(titleStr)
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        legend('pts ref','pts move');
        warning('on','MATLAB:legend:IgnoringExtraEntries')
        ax(axNum).Visible='on'; ax(axNum).Color='none';
        

    end
    
end




end