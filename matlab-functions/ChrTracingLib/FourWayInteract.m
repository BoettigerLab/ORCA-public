function  [probChange,numObs,typeM,distAcB,distAnB] = FourWayInteract(distMap,varargin)
% adapted from FourWayInteraction
% added conditioning against collapse for all 4 way pairs, not just 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'threshold','nonnegative',200};
defaults(end+1,:) = {'showPlots','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);


%%
[nReads,~,nCells] = size(distMap);
pts = 1:nReads;
npts = length(pts);
contactMap = distMap < pars.threshold; 

probChange = zeros(npts,npts,npts,npts);
numObs = zeros(npts,npts,npts,npts);
distAcB =  zeros(npts,npts,npts,npts);
distAnB =  zeros(npts,npts,npts,npts);
typeM = numObs;
tic
parfor ai=1:npts  % can parfor
    for bi=1:npts
        for ci=1:npts
            for di=1:npts
                a = pts(ai); 
                b = pts(bi);
                c = pts(ci);
                d = pts(di); 
                isCross1 = false;
                isCross2 = false;
                isInside = false;
                isOutside = false;
                isThreeWay = false;
                isSeparate = false;
                if (a<b && c<d) && ~(a == c && b== d)
                    if (c<a && d<b && d>a ) 
                        isCross1 = true;
                        typeM(ai,bi,ci,di) = 1;
                    elseif (c>a && c<b && d>b ) 
                        isCross2 = true;
                        typeM(ai,bi,ci,di) = 1;
                    elseif c>a && d<b  % c-d is inside a-b
                        isInside = true;
                        typeM(ai,bi,ci,di) = 2;
                    elseif c<a && d>b
                        isOutside = true;
                        typeM(ai,bi,ci,di) = 3;
                    elseif c==a || d==b || c==b %  || d==a
                        isThreeWay =true;
                        typeM(ai,bi,ci,di) = 4;
                    else
                        isSeparate = true; 
                        typeM(ai,bi,ci,di) = 5;
                    end
                         
                    hasData = ~isnan(distMap(a,b,:)) & ~isnan(distMap(c,d,:));
                    given = squeeze(contactMap(a,b,hasData)); 
                    condition = squeeze(contactMap(c,d,hasData));                    
                    extra =[];
                    if isCross1
                        extra =  squeeze(contactMap(a,c,hasData)) ;
                    elseif isCross2
                        extra =  squeeze(contactMap(b,d,hasData));
                    elseif isSeparate && c<a
                        extra =  squeeze(contactMap(a,c,hasData)) ;
                    elseif isSeparate && c>b    
                        extra =  squeeze(contactMap(b,d,hasData)) ;    
                    elseif isThreeWay 
                        if  c==a  
                            extra = squeeze(contactMap(c,d,hasData)) ;
                        elseif  d==b 
                            extra = squeeze(contactMap(d,c,hasData)) ;
                        elseif c==b
                            extra = squeeze(contactMap(b,d,hasData)) ;
                        end
                    elseif isInside || isOutside
                        extra = squeeze(contactMap(a,c,hasData)) ;
                    else
                        extra =[];
                    end
                 
                    if ~isempty(extra)
                        probChange(ai,bi,ci,di) =  sum(given & condition & ~extra)/sum(hasData) ./ ((sum(given)/sum(hasData))*(sum(condition)/sum(hasData))*(sum(~extra)/sum(hasData)) );
                        numObs(ai,bi,ci,di) = sum(given & condition  & ~extra);
                        distAcB(ai,bi,ci,di) = nanmedian(distMap(c,d,given & ~extra),3);
                        distAnB(ai,bi,ci,di) = nanmedian(distMap(c,d,~given & ~extra),3);
                    else
                        probChange(ai,bi,ci,di) =  sum(given & condition)/sum(hasData) ./ ((sum(given)/sum(hasData))*(sum(condition)/sum(hasData)));
                        numObs(ai,bi,ci,di) = sum(given & condition );
                        distAcB(ai,bi,ci,di) = nanmedian(distMap(c,d,given),3);
                        distAnB(ai,bi,ci,di) = nanmedian(distMap(c,d,~given),3);
                    end
                    
                end
            
            end
        end
    end
end
