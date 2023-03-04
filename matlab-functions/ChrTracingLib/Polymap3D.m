function [tform3D,xyzFix,fixDist] = Polymap3D(xyzRef,xyzWarp,varargin)
% [tform3D,xyzFix] = Polymap3D(xyzRef,xyzWarp,varargin)
% takes a two matched lists: xyzRef and xyzWarp. Both are Nx3 arrays of 
% x,y,z coordinates. The lists are ordered such that the ith element of
% each list corresponds to the same spot. The lists must be the same
% length.
% Returns tform3D containing the polynomial map,such that:
% xyzFix = tforminv(tform3D,xyzWarp(:,1),xyzWarp(:,2),xyzWarp(:,3));
% Also returns the result of this transform as an optional second argument
% For speed, if no plots are requested and only one output is requested,
% xyzFix will not be computed
% 
% -------------------------------
% Alistair Boettiger, Jan 2019 CC BY


defaults = cell(0,3);
defaults(end+1,:) = {'chnColors','colormap',hsv(2)};
defaults(end+1,:) = {'polyOrder','integer',2};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'figDist','integer',23};
defaults(end+1,:) = {'figMap','integer',22};
defaults(end+1,:) = {'max2D','nonnegative',0};
defaults(end+1,:) = {'max3D','nonnegative',0};
defaults(end+1,:) = {'maxPts','nonnegative',1e4};
defaults(end+1,:) = {'method',{'polynomial','affine'},'polynomial'};
defaults(end+1,:) = {'bins','array', 0:10:500};
defaults(end+1,:) = {'units','string', 'nm'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% max sep
if pars.max2D > 0
    origDist2D = sqrt(sum( (xyzWarp(:,1:2) - xyzRef(:,1:2)).^2,2));
    cut = origDist2D > pars.max2D;
    xyzWarp(cut,:) = [];
    xyzRef(cut,:) = [];
end
if pars.max3D > 0
    origDist3D = sqrt(sum( (xyzWarp(:,1:3) - xyzRef(:,1:3)).^2,2));
    cut = origDist3D > pars.max3D;
    xyzWarp(cut,:) = [];
    xyzRef(cut,:) = [];
end

%================== COMPUTE Shift ================%    
if strcmp(pars.method,'polynomial')
    tform3D = cp2tform3D(xyzRef,xyzWarp,'polynomial',pars.polyOrder);    
elseif strcmp(pars.method,'affine')
    tform3D = cp2tform3D(xyzRef,xyzWarp,'affine');    
end
if nargout > 1 || pars.showPlots
    xyzFix = tforminv(tform3D,xyzWarp(:,1),xyzWarp(:,2),xyzWarp(:,3));
    fixDist =sqrt(sum( (xyzWarp(:,1:3) - xyzRef(:,1:3)).^2,2));
else
    xyzFix = [];  % not defined if not plotting
    fixDist = [];
end


% ----- Just plotting
c=1;%  (there are only 2 colors);
if pars.showPlots
    nPts = size(xyzRef,1); % for speed
    if nPts > pars.maxPts
        s = randperm(nPts,pars.maxPts);
    else
        s = 1:nPts;
    end
    figure(pars.figMap); clf;
    subplot(2,2,1);
    x1 = xyzWarp(s,1); 
    y1 = xyzWarp(s,2);
    u1 = xyzFix(s,1)-xyzWarp(s,1); 
    v1 = xyzFix(s,2)-xyzWarp(s,2); 
    plot(xyzWarp(s,1),xyzWarp(s,2),'o','color',pars.chnColors(c+1,:)); hold on;
    plot(xyzFix(s,1),xyzFix(s,2),'.','color',pars.chnColors(c+1,:),'MarkerSize',15);
    plot(xyzRef(s,1),xyzRef(s,2),'+','color',pars.chnColors(1,:)); hold on;
    % add connecting lines
    plot([xyzRef(s,1),xyzWarp(s,1)]',[xyzRef(s,2),xyzWarp(s,2)]','-','color',pars.chnColors(1,:)); hold on;
    plot([xyzRef(s,1),xyzFix(s,1)]',[xyzRef(s,2),xyzFix(s,2)]','-','color',pars.chnColors(1,:)); hold on;
    legend('orig','corrected','target');
    title('xy projection');
    subplot(2,2,2);
    x2 = xyzWarp(s,1); 
    y2 = xyzWarp(s,3);
    u2 = (xyzFix(s,1)-xyzWarp(s,1)); 
    v2 = (xyzFix(s,3)-xyzWarp(s,3)); 
    plot(xyzWarp(s,1),xyzWarp(s,3),'o','color',pars.chnColors(c+1,:)); hold on;
    plot(xyzFix(s,1),xyzFix(s,3),'.','color',pars.chnColors(c+1,:),'MarkerSize',15);
    plot(xyzRef(s,1),xyzRef(s,3),'+','color',pars.chnColors(1,:)); hold on;
     % add connecting lines
    plot([xyzRef(s,1),xyzWarp(s,1)]',[xyzRef(s,3),xyzWarp(s,3)]','-','color',pars.chnColors(1,:)); hold on;
    plot([xyzRef(s,1),xyzFix(s,1)]',[xyzRef(s,3),xyzFix(s,3)]','-','color',pars.chnColors(1,:)); hold on;
    legend('orig','corrected','target');
    title('xz projection');
    subplot(2,2,3);
    quiver(x1,y1,u1,v1); 
    title('xy projection, vector field');
    subplot(2,2,4);
    quiver(x2,y2,u2,v2*.04); hold on;
     title('xz projection, vector field');
    set(gcf,'color','w');
    

    bins = pars.bins; % 
    stp = bins(2)-bins(1);
    origDist = sqrt(sum( (xyzWarp(:,1:2) - xyzRef(:,1:2)).^2,2));
    fixDist = sqrt(sum( (xyzFix(:,1:2) - xyzRef(:,1:2)).^2,2));
    figure(pars.figDist); clf; 
    subplot(3,2,2); hist(origDist,bins);  xlim([0,max(bins)+stp]); title(['Orig XY dist e=',num2str(median(origDist),4)]); xlabel(pars.units);
    subplot(3,2,4); hist(fixDist,bins);  xlim([0,max(bins)+stp]); title(['Cor. XY dist e=',num2str(median(fixDist),4)]); xlabel(pars.units);
    origDist = sqrt(sum( (xyzWarp(:,1:3) - xyzRef(:,1:3)).^2,2));
    fixDist = sqrt(sum( (xyzFix(:,1:3) - xyzRef(:,1:3)).^2,2));
    subplot(3,2,1); hist(origDist,bins); xlim([0,max(bins)+stp]); title(['Orig 3D dist e=',num2str(median(origDist),4)]); xlabel(pars.units);
    subplot(3,2,3); hist(fixDist,bins);  xlim([0,max(bins)+stp]);  title(['Cor. 3D dist e=',num2str(median(fixDist),4)]); xlabel(pars.units);
    subplot(3,2,5); hist(abs(u1)+abs(v1),bins);  xlim([0,max(bins)+stp]); title('xy shift mag'); xlabel(pars.units);
    subplot(3,2,6); hist(abs(v2),bins);  xlim([0,max(bins)+stp]); title('z shift mag'); xlabel(pars.units);
    set(gcf,'color','w');
end