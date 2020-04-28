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
defaults(end+1,:) = {'max2D','nonnegative',0};
defaults(end+1,:) = {'max3D','nonnegative',0};
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
tform3D = cp2tform3D(xyzRef,xyzWarp,'polynomial',pars.polyOrder);    
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
    figure(22); clf;
    subplot(1,2,1);
    x = xyzWarp(:,1); 
    y = xyzWarp(:,2);
    u = xyzFix(:,1)-xyzWarp(:,1); 
    v = xyzFix(:,2)-xyzWarp(:,2); 
    quiver(x,y,u,v); hold on;
    plot(xyzWarp(:,1),xyzWarp(:,2),'o','color',pars.chnColors(c+1,:)); hold on;
    plot(xyzFix(:,1),xyzFix(:,2),'.','color',pars.chnColors(c+1,:),'MarkerSize',15);
    plot(xyzRef(:,1),xyzRef(:,2),'+','color',pars.chnColors(1,:)); hold on;
    legend('quiver','orig','corrected','target');
    title('xy projection');
    subplot(1,2,2);
    x = xyzWarp(:,1); 
    y = xyzWarp(:,3);
    u = (xyzFix(:,1)-xyzWarp(:,1)); 
    v = (xyzFix(:,3)-xyzWarp(:,3)); 
    quiver(x,y,u,v,0); hold on;
    plot(xyzWarp(:,1),xyzWarp(:,3),'o','color',pars.chnColors(c+1,:)); hold on;
    plot(xyzFix(:,1),xyzFix(:,3),'.','color',pars.chnColors(c+1,:),'MarkerSize',15);
    plot(xyzRef(:,1),xyzRef(:,3),'+','color',pars.chnColors(1,:)); hold on;
    legend('quiver','orig','corrected','target');
    title('xz projection');
    

    bins = 0:10:500;
    origDist = sqrt(sum( (xyzWarp(:,1:2) - xyzRef(:,1:2)).^2,2));
    fixDist = sqrt(sum( (xyzFix(:,1:2) - xyzRef(:,1:2)).^2,2));
    figure(23); clf; 
    subplot(2,2,2); hist(origDist,bins);  xlim([0,max(bins)+1]); title(['Orig XY dist e=',num2str(median(origDist),4)]);
    subplot(2,2,4); hist(fixDist,bins);  xlim([0,max(bins)+1]); title(['Cor. XY dist e=',num2str(median(fixDist),4)]);
    xlabel('nm');
    origDist = sqrt(sum( (xyzWarp(:,1:3) - xyzRef(:,1:3)).^2,2));
    fixDist = sqrt(sum( (xyzFix(:,1:3) - xyzRef(:,1:3)).^2,2));
    subplot(2,2,1); hist(origDist,bins); xlim([0,max(bins)+1]); title(['Orig 3D dist e=',num2str(median(origDist),4)]);
    subplot(2,2,3); hist(fixDist,bins);  xlim([0,max(bins)+1]);  title(['Cor. 3D dist e=',num2str(median(fixDist),4)]);
    xlabel('nm');
end