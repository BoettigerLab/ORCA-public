function [insScore,sepScore,sepRatio] = PolymerSepPlane(polymerAll,varargin)
% [insScore,sepScore] = PolymerSepPlane(polymerAll,varargin)
% insScore is number of cross plane localizations / total localizations
% sepScore is ave local distance (box on diagonal) over ave distal distance
%   (box shifted one box length off the diagonal).  

defaults = cell(0,3);
defaults(end+1,:) = {'w','integer',5};
defaults(end+1,:) = {'parallel','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'distMap','array',[]};
defaults(end+1,:) = {'computeSepPlane','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

w = pars.w;
verbose = pars.verbose;
computeSepPlane = pars.computeSepPlane;
% polymerAll = cat(3,polymerFilts{:});

nCells = size(polymerAll,3); 
nReads = size(polymerAll,1); 
insScore = nan(nReads,nCells);
sepScore = nan(nReads,nCells);
sepRatio = nan(nReads,nCells);
distMap = pars.distMap; 
tic
if pars.parallel
    parfor k=1:nCells
        string = polymerAll(:,1:3,k); % 4, 12
        % figure(1); clf; 
        % PlotPolymerTube(string);
        for r=1:nReads  % r = 33
           % if r>=w && r<=nReads-w+1
                s1 = max(1,r-w);
                s2 = min(nReads,r+w);
                setA = string(s1:r-1,:);
                setB = string(r:s2,:);

                setA( isnan(setA(:,1)) ,:) = [];
                setB( isnan(setB(:,1)) ,:) = [];               
                if size(setA,1) > 1 && size(setB,1) > 1 && computeSepPlane
                    try
                        [cA,cB] = SepPlane(setA,setB,'showPlots',false);
                        % nansum([cA;cB])/length([cA; cB])
                        insScore(r,k) = nansum([cA;cB])/length([cA; cB]);
                    catch er
                        disp('problem');
                    end    
                end
                if ~isempty(distMap)
                    intReg = distMap(s1:r-1,r:s2,k);
                    intReg(intReg==0) = NaN;
                    sepScore(r,k) = nanmean(intReg(:));
                    
                    near = distMap(s1:s2,s1:s2,k);
                    near(near==0) = NaN;
                    sM = min(nReads,s2+2*w);
                    far1 = distMap(s1:s2,s2+1:sM,k);
                    sM = max(1,s1-2*w);
                    far2 = distMap(s1:s2,sM:s1-1,k);
                    far = [far1(:); far2(:)];
                    far(far==0) = NaN;
                    sepRatio(r,k) = nanmedian(far(:))/nanmedian(near(:));  
                end
            % end
        end
        if verbose
            disp(['processing cell ',num2str(k)]);
        end
    end
    
else
    for k=1:nCells
        string = polymerAll(:,1:3,k); % 4, 12
        % figure(1); clf; 
        % PlotPolymerTube(string);
        for r=1:nReads  % r = 33
            if r>=w && r<=nReads-w+1
                s1 = max(1,r-w);
                s2 = min(nReads,r+w);
                setA = string(s1:r-1,:);
                setB = string(r:s2,:);

                setA( isnan(setA(:,1)) ,:) = [];
                setB( isnan(setB(:,1)) ,:) = [];               
                if size(setA,1) > 1 && size(setB,1) > 1 && computeSepPlane
                    try
                        [cA,cB] = SepPlane(setA,setB,'showPlots',false);
                        % nansum([cA;cB])/length([cA; cB])
                        insScore(r,k) = nansum([cA;cB])/length([cA; cB]);
                    catch er
                        disp('problem');
                    end    
                end
                if ~isempty(distMap)
                    intReg = distMap(s1:r-1,r:s2,k);
                    intReg(intReg==0) = NaN;
                    sepScore(r,k) = nanmean(intReg(:));

                    near = distMap(s1:s2,s1:s2,k);
                    near(near==0) = NaN;
                    sM = min(nReads,s2+2*w);
                    far1 = distMap(s1:s2,s2+1:sM,k);
                    sM = max(1,s1-2*w);
                    far2 = distMap(s1:s2,sM:s1-1,k);
                    far = [far1(:); far2(:)];
                    far(far==0) = NaN;
                    sepRatio(r,k) = nanmedian(far(:))/nanmedian(near(:));  
                end
            end
        end
        if verbose
            disp(['processing cell ',num2str(k)]);
        end
    end
    
end

if verbose
    t = toc;
    disp(['processing completed in ' num2str(t/60), ' min']);
end

