function NPolyNeighbors(distMaps,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'badHybes','integer',[]};
defaults(end+1,:) = {'threshold','positive',150};
pars = ParseVariableArguments(varargin,defaults,mfilename);

theta = pars.threshold;
[nB,~,nCells] = size(distMaps);
badHybe = false(1,nB);
badHybe(pars.badHybes) = true;
neib = zeros(nB,nCells);
neibID = zeros(nB,nCells);
for c = 1:nCells% c =20
    map = distMaps(:,:,c);
      %   figure(1); clf; imagesc(map); colorbar; caxis([0,1200]); pause(.1);
    for b=1:nB % b = 4
        neibs = map(b,:) < theta;
        neib(b,c) = sum(neibs) ; % ./ sum(map(b,:) < inf);
        neibID(b,c) = median( find( neibs) )-b;
    end
    skip = isnan(map(1,:)) | badHybe;
    neib(skip,c) = nan;
    neibID(skip,c) = nan;
end
% covMap = sum(~isnan(distMaps),3)./nCells; figure(4); clf; imagesc(covMap);
% cov = nanmean(covMap,2);
m = nanmean( neib,2);
figure(2); hold on; 
sem = nanstd(neib,[],2)./sqrt( sum(~isnan(neib),2) );
% plot(m+2*sem,'k--'); hold on;
% plot(m-2*sem,'k--');
plot(m,'-','LineWidth',2);


figure(3); hold on;
plot( nanmean(neibID,2));