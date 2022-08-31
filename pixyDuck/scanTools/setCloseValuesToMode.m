function allRoundRes=setCloseValuesToMode(allRoundRes)
% can have issue where resolution is eps(res) away from another resolution value. Just pick most popular one (median) and use this.

assert(isvector(allRoundRes))
if size(allRoundRes,2)==1
    allRoundRes=allRoundRes'; % make row vector
end

%allRoundRes=[0.6+100*eps(0.6),0.1,allRoundRes,0.6,0.7,0.1+100*eps(1)]
%numUniqueIn=length(unique(allRoundRes))

epsMultiplier=1000; % willing to tolerate up to epsMultiplier*eps(max(Resolutions) distance away

[uniqueResolutionsRaw,indFirstUniqueResRawAll,indRoundOfRes]=unique(allRoundRes);

for indThis=indFirstUniqueResRawAll'
    absDiffVect=abs(allRoundRes-allRoundRes(indThis))';
    indTinyDiffFromThisRes=find(all([absDiffVect>0,absDiffVect<(epsMultiplier*eps(max(uniqueResolutionsRaw)))],2));
    if ~any(indTinyDiffFromThisRes)
        % move on to next iResRaw
    else
        % set all of these rounds to be equal to the most popular one
        indSelfAndTinyDiffOthers=[indThis;indTinyDiffFromThisRes];
        % get mode of these
        resMode=mode(allRoundRes(indSelfAndTinyDiffOthers));
        % upated current list of allRoundRes to this mode
        allRoundRes(indSelfAndTinyDiffOthers)=resMode;
        % we will always return to the next 'unique'(but not really)
        % resolutions in the for loop but it won't find any difference
        % now and won't reset it
    end
end


%allRoundRes
%numUniqueOut=length(unique(allRoundRes))





end