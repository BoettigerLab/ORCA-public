function [probChange,numObs,typeM] = MultiLoopEffects(distMap,theta)
%% updates
% 
%% To Do
% - let's try a version where we remove the 'collapsed state' (a,b,c,d
% <theta) from all data, rather than just conditioning against it in the
% cross-loops, (which likely over-estimates the exclusion of cross-loops). 

%% loop interactions

npts = size(distMap,1);
pts = 1:npts;
contactMap = distMap < theta; 
[nReads,~,nCells] = size(contactMap);

probChange = zeros(npts,npts,npts,npts);
numObs = zeros(npts,npts,npts,npts);
typeM = numObs;
tic
parfor ai=1:npts  % can parfor
    for bi=1:npts % ai:npts (doesnt' work in parfor
        for ci=1:npts
            for di=1:npts % ci:npts
                a = pts(ai); 
                b = pts(bi);
                c = pts(ci);
                d = pts(di); 
                isCross1 = false;
                isCross2 = false;
                isThreeWay = false;
                if (a<b && c<d) && ~(a == c && b== d) % this enforces the upper diagonal
                    if (c<a && d<b && d>a ) %   d  
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
                    elseif c==a || d==b || c==b || d==a 
                        isThreeWay =true;
                        typeM(ai,bi,ci,di) = 4;
                    else
                        isSeparate = true; 
                        typeM(ai,bi,ci,di) = 5;
                    end
                         
                    hasData = ~isnan(distMap(a,b,:)) & ~isnan(distMap(c,d,:));
                    given = squeeze(contactMap(a,b,hasData)); 
                    condition = squeeze(contactMap(c,d,hasData));                    

                    if isCross1
                        extra = squeeze(contactMap(a,c,hasData)) ;
                    elseif isCross2
                        extra = squeeze(contactMap(b,d,hasData));
                    else
                        extra =[];
                    end
                 
                    % if ab what happens to cd
                    %   for 3-way  compare ba bc to ac
                    %     or   P(ab | bc)/P(ab)P(bc) 
                    %    NOT   P(ac | ab) for a<c<b
                    if isThreeWay
                        bc = b==c && a < b && c < d;
                        % ac = a==c && b < a && a < d;
                        % ad = a==d && b < a && a < c; % ba ac
                        % bd = b==d && a < b && d < c;
                        if  bc %  || ac || ad || bd
                            % if a==3 && b==10
                            %     disp('found 3-way')
                            % end
                            probChange(ai,bi,ci,di) =  sum(given & condition)/sum(hasData) ./ ((sum(given)/sum(hasData))*(sum(condition)/sum(hasData)));
                            numObs(ai,bi,ci,di) = sum(given & condition );
                        end
                    else
                        if ~isempty(extra)
                            probChange(ai,bi,ci,di) =  sum(given & condition & ~extra)/sum(hasData) ./ ((sum(given)/sum(hasData))*(sum(condition)/sum(hasData))*(sum(~extra)/sum(hasData)) );
                            numObs(ai,bi,ci,di) = sum(given & condition  & ~extra);
                        else
                            probChange(ai,bi,ci,di) =  sum(given & condition)/sum(hasData) ./ ((sum(given)/sum(hasData))*(sum(condition)/sum(hasData)));
                            numObs(ai,bi,ci,di) = sum(given & condition );
                        end
                    end
                    
                end
            
            end
        end
    end
end

toc
%%
% % a = 4; b=9; %  very good ex for IMR90
% a = 3; b=10; 
% mapN = squeeze(numObs(a,b,:,:));
% 
% map = log2( squeeze(probChange(a,b,:,:)) );
% map(mapN<6) = 0;
% probMapFig = figure(3); clf;
% imagesc(map); colorbar;
% GetColorMap('RedWhiteBlue');
% caxis([-2,2]);
% hold on; plot(b,a,'m.','MarkerSize',15); axis image;
% typeFig = figure(4); clf; imagesc( squeeze(typeM(a,b,:,:))); 
% colormap(magma(6)); colorbar; axis image;