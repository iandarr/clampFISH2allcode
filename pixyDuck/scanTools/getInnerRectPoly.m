function poly=getInnerRectPoly(polyvec,varargin)
assert(isa(polyvec,'polyshape'))

if nargin==2
    if islogical(varargin{1})
        plotOn=varargin{1};
    end
elseif nargin>2
    error('too many inputs')
else
    plotOn=false;
end
    
polyIn=union(polyvec);
polyIn=polyIn.rmholes;


poly=polyIn; % poly is working polygon in loop

[tempx,tempy]=poly.boundary; % pts are clockwise
xy=[tempx(1:end-1) tempy(1:end-1)]; % last point is same as first point, don't count it.
clear tempx tempy
np=size(xy,1);
npStart=np;

if plotOn
        fh=figure; clf;
        ax=gca; hold on;    
    plot(polyIn,'FaceColor','b','FaceAlpha',1)
	labels=arrayfun(@num2str,1:np-1,'UniformOutput',0)';
	hlabels=text(xy(1:end-1,1),xy(1:end-1,2),labels)
    %text(xy(end,1)+diff(ax.XLim)/30,xy(end,2)+diff(ax.YLim)/30,num2str(size(xy,1)))

end


    
finalPts=4;

temp=[polyvec.Vertices];
xyOrigPolyvec=[reshape(temp(:,1:2:end),4*length(polyvec),1),reshape(temp(:,2:2:end),4*length(polyvec),1)];


while np>finalPts % np (num points) = num sides
    
    
    %% find sides. poly goes around clockwise direction
    sp1=xy(1:end,:);
    sp2=[xy(2:end,:);xy(1,:)];
    % figure out if side is horizontal or vertical for later
    isHorizontal=sp1(:,2)==sp2(:,2);
    isVertical=sp1(:,1)==sp2(:,1);
    
    if ~all((isHorizontal+isVertical)==1) % are all sides vertical and horizontal?
        
        % first try cleaning data, since
        % sometimes after union, sides that should be horizontal/vertical
        % have a rise/run within eps(sp1) or so

        idxNotHorV=~((isHorizontal+isVertical)==1);
        indNotHorV=find(idxNotHorV);
        
        xy=cleanXY(xy,xyOrigPolyvec);
        
        % update sp1 and sp2
        sp1=xy(1:end,:);
        sp2=[xy(2:end,:);xy(1,:)];
        % figure out if side is horizontal or vertical for later
        isHorizontal=sp1(:,2)==sp2(:,2);
        isVertical=sp1(:,1)==sp2(:,1);
        
        if ~all((isHorizontal+isVertical)==1)
            temp=join(cellstr(num2str(indNotHorV))',', ');
            error('there are some sides on the perimeter that are not horizontal or vertical. Must be for this function to work. they are %s',temp{1});
        end
        
    end
    %% if we extend each side in both directions and try to make a new smaller polygon, which points are now outside the new polygon
    % right-side normal vect for each side. Not a unit vector
    %cameraAngle=-90;
    %R = [cosd(cameraAngle) -sind(cameraAngle); sind(cameraAngle) cosd(cameraAngle)];
    %R=[ 0 1; -1 0]; % for right-side normal vect
    R=[ 0 -1; 1 0]; % for right-side normal vect
    Ns=[R*[sp2-sp1]']';
    NsStack=reshape(repmat(Ns',np,1),2,np*np)';
    sp2ToPtsX=(xy(:,1)-repmat(sp1(:,1)',np,1))';
    sp2ToPtsY=(xy(:,2)-repmat(sp1(:,2)',np,1))';
    sp2ToPtsXstack=reshape(sp2ToPtsX,1,np*np)';
    sp2ToPtsYstack=reshape(sp2ToPtsY,1,np*np)';
    sp2ToPtsstack=[sp2ToPtsXstack,sp2ToPtsYstack];
    vToRightStack=dot(NsStack',sp2ToPtsstack')'; % doesn't matter whether we subtract s1 or s2 since it's normal
    
    
    sOutStack=vToRightStack<0; % numRem = num to remove
    sOut=reshape(sOutStack,[np,np])';
    
    
    
    %% FORWARD direction
    % extending sides forward
    sOutNext3Fwd=sOutNextN(sOut,3,'forward');
    idxsCandFwd=all(sOutNext3Fwd==[1 1 0],2);
    indsCandFwd=find(idxsCandFwd);
    numCandFwd=length(indsCandFwd);
    % for each candidate side to extend, get the rectangle that would be removed so we can find the smallest of these
    % rectangles' areas, and that side will be chosen for extension
    %
    % we know that for a side index s (with points s and s+1), when extended, it intersects between points s+3 and s+4 points
    %xyNewFwd=polyxpoly
    
    % using vectors with height of candidates
    isH=isHorizontal(idxsCandFwd);
    isV=isVertical(idxsCandFwd);
    
    sp1CandFwd=sp1(idxsCandFwd,:);
    
    xNewFwd=nan(numCandFwd,1);
    xNewFwd(isH)=xy(wrap(indsCandFwd(isH)+3,np),1);
    xNewFwd(isV)=sp1CandFwd(isV,1);
    yNewFwd=nan(numCandFwd,1);
    yNewFwd(isH)=sp1CandFwd(isH,2);
    yNewFwd(isV)=xy(wrap(indsCandFwd(isV)+3,np),2);
    
    xyNewFwd=[xNewFwd,yNewFwd];
    
    % xy of areas to be removed (forward)
    xRemFwd=[xy(wrap(indsCandFwd+1,np),1),xy(wrap(indsCandFwd+2,np),1),xy(wrap(indsCandFwd+3,np),1),xyNewFwd(:,1)];
    yRemFwd=[xy(wrap(indsCandFwd+1,np),2),xy(wrap(indsCandFwd+2,np),2),xy(wrap(indsCandFwd+3,np),2),xyNewFwd(:,2)];
    
    %% BACKWARDS direction
    % extending sides backwards (sOut 3 columns also is left-->right looking towards backwards)
    sOutNext3Bck=sOutNextN(sOut,3,'backwards');
    idxsCandBck=all(sOutNext3Bck==[1 1 0],2);
    indsCandBck=find(idxsCandBck);
    numCandBck=length(indsCandBck);
    
    % using vectors with height of candidates
    isH=isHorizontal(idxsCandBck);
    isV=isVertical(idxsCandBck);
    
    sp1CandBck=sp1(idxsCandBck,:);
    
    xNewBck=nan(numCandBck,1);
    xNewBck(isH)=xy(wrap(indsCandBck(isH)-2,np),1);
    xNewBck(isV)=sp1CandBck(isV,1);
    yNewBck=nan(numCandBck,1);
    yNewBck(isH)=sp1CandBck(isH,2);
    yNewBck(isV)=xy(wrap(indsCandBck(isV)-2,np),2);
    
    xyNewBck=[xNewBck,yNewBck];
    
    % xy of areas to be removed (forward)
    xRemBck=[xyNewBck(:,1),xy(wrap(indsCandBck-2,np),1),xy(wrap(indsCandBck-1,np),1),xy(indsCandBck,1)];
    yRemBck=[xyNewBck(:,2),xy(wrap(indsCandBck-2,np),2),xy(wrap(indsCandBck-1,np),2),xy(indsCandBck,2)];
    
    
    %% remove smallest area polygon
    polyCandRem=polyshapeVector([xRemFwd;xRemBck],[yRemFwd;yRemBck])';
    polyArea=[polyCandRem.area];
    [~,indSmallest]=min(polyArea);
    polyRem=polyCandRem(indSmallest);
    
    polyNew=subtract(poly,polyRem);
    
     %% update poly, xy, and np
    poly=polyNew;
    [tempx,tempy]=poly.boundary; % pts are clockwise
    xy=[tempx(1:end-1) tempy(1:end-1)]; % last point is same as first point, don't count it.
    clear tempx tempy
    np=size(xy,1);
    
    if np<=finalPts
        xy=cleanXY(xy,xyOrigPolyvec); % may need to clean
        poly=polyshape(xy(:,1)',xy(:,2)');
    end
    
    %% plot if that's on
    if plotOn
        hrem=plot(polyRem,'FaceAlpha',1,'FaceColor','r');
        
        
        if 5000 * polyRem.area < abs(diff(ax.XLim)*diff(ax.YLim))
           % very small, create a brief glow around it
           glowDist=min(abs(diff(ax.XLim)),abs(diff(ax.YLim)))/100;
           polyRemGlow=polybuffer(polyRem,glowDist);
           hGlow=plot(polyRemGlow,'FaceAlpha',0.2,'FaceColor','r');
           pause(1/sqrt(npStart)) % to show main face
           delete(hGlow)
        end
        pause(2/sqrt(npStart)) % to show main face
        
        delete(hrem)
        if exist('hpoly','var')
            delete(hpoly)
        end
        hpoly=plot(polyNew,'FaceAlpha',0.8,'FaceColor','g');
        
        labels=arrayfun(@num2str,1:np,'UniformOutput',0)';
        if exist('hlabels','var')
            delete(hlabels)
        end
        hlabels=text(xy(1:end,1),xy(1:end,2),labels);

    end
    

    
end


end

function indOut=wrap(ind,n)
% wrap indices to vector of length n with rollover back to one if you reach
% the end.
indOut= mod(ind-1, n) + 1;
end




function sOutNextN=sOutNextN(sOut,n,direction)
assert(ismember(direction,{'forward','backwards'}))

assert(n<=size(sOut,2))
np=size(sOut,1);
sOutNextN=nan(np,1);

switch direction
    case 'forward'
        sOutPad=[sOut,sOut(:,1:n+1)];
        for i=1:np
            sOutNextN(i,1:n)=sOutPad(i,(i+2):(i+n+1));
            %numRemNext3=arrayfun(@(x) x(x:x+n-1),2),[sRem,sRem(:,1:n)])
        end
        
    otherwise %'backwards'
        sOutPad=[sOut(:,np-n+1:np),sOut];
        for i=1:np
            sOutNextN(i,1:n)=sOutPad(i,(n+i-1):-1:(i));
            %numRemNext3=arrayfun(@(x) x(x:x+n-1),2),[sRem,sRem(:,1:n)])
        end
end
end


function xy=cleanXY(xy,xyOrig)
   np=size(xy,1);
        sp1=xy(1:end,:);
            sp2=[xy(2:end,:);xy(1,:)];
isHorizontal=sp1(:,2)==sp2(:,2);
isVertical=sp1(:,1)==sp2(:,1);
idxNotHorV=~((isHorizontal+isVertical)==1);
indNotHorV=find(idxNotHorV);

sideDirs=sp2(idxNotHorV,:)-sp1(idxNotHorV,:);
        for i=1:length(indNotHorV)
            ind=indNotHorV(i);
            sideDir=sideDirs(i,:);
            
            [absErr,problemDim]=min(abs(sideDir(sideDir~=0)));
            epsAtPt=eps(xy(ind,problemDim));
            tolerance= 200 * epsAtPt;
            if absErr<tolerance
                % change value to closest original point value if it's
                % within tolerance
                replaceIdx=abs(xy(ind,problemDim)- xyOrig(:,problemDim)) < tolerance;
                indXYreplace=find(replaceIdx,1); % use first one
                if ~isempty(indXYreplace)
                    valReplace=xyOrig(indXYreplace,problemDim);
                    assert(abs(xy(ind,problemDim)-valReplace)<tolerance) % unecessary but checking my logic here
                    xy(ind,problemDim)=valReplace;                    
                    % although we don't know the other point in segment
                    % isn't valReplace, we might as well change it.
                    if ind<np
                        xy(ind+1,problemDim)=valReplace;
                    else % ind == np
                        xy(1,problemDim)=valReplace;
                    end
                else
                    % we couldn't find an original point dimension close enough to
                    % serve as a template.
                    error('hmmm why does this happen with rectangular data')
                end
                else
                % else do nothing, can't clean it
                warning('getInnerRectPoly: cannot clean a data point to get the inner rect')
            end
            
        end

end