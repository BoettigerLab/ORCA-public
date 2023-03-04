function [entropy_nt, bootstrap]= ChainEntropy(xyz,x,varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'badHybes','positive',[]};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'bootstrap','positive',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);


[nB,~,nC] = size(xyz);
entropy_nt = NaN(nB,nB,pars.bootstrap);
for b=1:pars.bootstrap
    if pars.bootstrap > 1
        resample = randi(nC,nC,1);
    else
        resample = 1:nC;
    end
    xyzR = xyz(:,:,resample);
    for ref =1:nB
        dXYZ = abs(xyzR(:,:,:) -repmat(xyzR(ref,:,:),nB,1,1) );
        if ~isempty(pars.badHybes)
            dXYZ(pars.badHybes-1,:,:) = NaN;
        end
        dD = squeeze(sqrt(sum(dXYZ.^2,2)));   % xyz coverted to dist 
        for r=1:size(dD,1) % 
            [p,x] = hist(dD(r,:),x); % compare varation across cells
            % h = histogram(dD(r,:),x);  % MUCH slower
%             p = h.BinCounts;
            if r==28 && pars.showPlots  && b==1
                figure(4); hold on; bar(x,p./sum(p));
            end
            p=p./nansum(p);
            entropy_nt(ref,r,b) = -nansum(p.*log2(p));
        end
        entropy_nt(ref,[ref,pars.badHybes-1],b) = nan;
        if pars.showPlots && b==1
            figure(2);   plot(entropy_nt(ref,:,b),'.-');
        end
        bootstrap = [];
    end

end
%% old version, diff nearest 
% 
% dXYZ = diff(xyz,1);
% if ~isempty(pars.badHybes)
% dXYZ(pars.badHybes-1,:,:) = NaN;
% end
% dD = squeeze(sqrt(sum(dXYZ.^2,2))); 
% nStp = size(dD,1);
% entropy_nt = NaN(nStp,1);
% for r=1:size(dD,1)
%     [p,x] = hist(dD(r,:),x);
%     if r==1 && pars.showPlots
%         figure(4); hold on; bar(x,p./sum(p));
%     end
%     p=p./nansum(p);
%     entropy_nt(r) = -nansum(p.*log2(p));
% end
% if pars.showPlots
%     figure(2);   plot(entropy_nt,'.-');
% end
% bootstrap = [];
% 
% if pars.bootstrap
%     [nStp,~,nCells] = size(xyz);
%     draws = randi(nCells,[nCells,pars.bootstrap]);
%     bootstrap = NaN(nStp,pars.bootstrap);
%     for b=1:pars.bootstrap
%         XYZ = xyz(:,:,draws(:,b));
%         dXYZ = diff(XYZ,1);
%         if ~isempty(pars.badHybes)
%             dXYZ(pars.badHybes-1,:,:) = NaN;
%         end
%         dD = squeeze(sqrt(sum(dXYZ.^2,2))); 
%         for r=1:size(dD,1)      
%             [p,x] = hist(dD(r,:),x);
%             p=p./nansum(p);
%             bootstrap(r,b) = -nansum(p.*log2(p));
%         end
%     end
% end
% 
% 
% entropy_nt(entropy_nt==0) = NaN;

