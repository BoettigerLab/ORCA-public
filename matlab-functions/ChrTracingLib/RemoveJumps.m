function [polymerOut,badPts] = RemoveJumps(polymerIn,varargin)
% Remove jumpy points from polymer path
%
%% Inputs
% polymerIn - required Nx3 polymer path  (Nx2 will also work)
% 
%% Optional inputs
% localRegion - number of upstream and downstream points to consider in determining if a jump is too big 
% maxDiff - max relative step size, recorded as fold change greater than the median step size of other points from their local area which is acceptable
% maxAbsStep - max absolute step size, in 
%% Outputs
% polymerOut - bad points replaced with NaNs
% badPts - (logical/boolean) indices of the coordinates to remove

defaults = cell(0,3);
defaults(end+1,:) = {'localRegion','integer',4}; % number of points in front and behind to use for estimating expected position 
defaults(end+1,:) = {'maxDiff','positive',2}; % fold change greater than the median step size of other points from their local area which is acceptable
defaults(end+1,:) = {'maxAbsStep','positive',inf}; % maximum absolute step size
defaults(end+1,:) = {'removeLoners','boolean',false}; % if a point is the only one within its local region, drop it. 
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

polymerOut = polymerIn;
nB = size(polymerOut,1);
if nB < 100
    m1s = squareform(pdist(polymerOut));  % too big for large polymers
    m1s = m1s - triu(m1s,pars.localRegion)  - tril(m1s,-pars.localRegion);
else % looping is slow for small polymers but necessary and faster for lager ones 
    m1s = zeros(nB,nB);
    for b=1:nB-1
        b1 = b; % (b-1)*pars.localRegion+1;
        b2 = min(nB,b+pars.localRegion-1); % b*pars.localRegion;
        m1s(b1:b2,b1:b2) = squareform(pdist(polymerOut(b1:b2,:)));
    end
end
% figure(7); clf; imagesc(m1s); colorbar; caxis([0,1000])

if pars.removeLoners
    badPts = nansum(m1s) == 0; %#ok<NANSUM>
else
    badPts = false(1,size(polymerIn,1));
end

m1s(m1s==0) = nan;
% figure(7); clf; imagesc(m1s); colorbar; caxis([0,1000])
m1 = nanmean(m1s);
badPts = m1 > nanmedian(m1)*pars.maxDiff | m1 > pars.maxAbsStep | badPts;
if pars.verbose
    disp(['removed points ',num2str(find(badPts))]);
end
 % figure(7); clf;  bar(m1);
polymerOut(badPts,:,1) = nan;

