function [pks,pkX,logRs,logNmeans,diffLogNmeans] = FractalDim(X)

% % % Simulated data to produce X
% meandist = 25;
% Ncluster = 500;
% X = arrayfun(@(x)[randn(1,Ncluster)+rand*meandist;randn(1,Ncluster)+rand*meandist;randn(1,Ncluster)+rand*meandist],1:20,'UniformOutput',false);
% X=cell2mat(X)';
% figure(1); clf; scatter3(X(:,1),X(:,2),X(:,3),'.'); PresentationPlot();

try
    
%      while length(X) > 3000
%          X = X(1:2:end,:);
%      end
    
% Actual Run
D = pdist(X);  % euclidean distance
D_sq = squareform(D);

%limits of distances
numPts = 30; 
logRs = linspace(log(min(D)),log(max(D)),numPts);
logNmeans = arrayfun(@(R)log(mean(sum(D_sq < exp(R),1) )),logRs);

dx = (logRs(end) - logRs(1) )/ numPts;
diffLogNmeans = diff(logNmeans)/dx;
[pks,locs] = findpeaks(diffLogNmeans,'MINPEAKHEIGHT',1);
pkX = logRs(locs);
catch er
    warning(er.getReport);
    pks = NaN;
    pkX = NaN;
    logRs = NaN;
    logNmeans= NaN;
    diffLogNmeans  = NaN;
end


% figure(4); clf; plot(logRs(1:end-1),diffLogNmeans,'b.-'); hold on;
% plot(pkX,pks,'k^');

% load sunspot.dat
% relNums=sunspot(:,2);
% figure(10); clf;
% plot(sunspot(:,2)); hold on;
% [pks,locs] = findpeaks(relNums);
% plot(locs,pks,'k.','MarkerSize',20);
% xlabel('Year');
% ylabel('Sunspot Number')
% title('Find All Peaks');
% 
% findpeaks(logRs',logNmeans')
% 
% figure(2); clf; plot(logRs,logNmeans)
% title('logNmeans vs logRs'); PresentationPlot();
% 
% figure(3); clf; plot(logRs(2:end),diff(logNmeans));
% title('d(logNmeans)/d(logRs)'); PresentationPlot();