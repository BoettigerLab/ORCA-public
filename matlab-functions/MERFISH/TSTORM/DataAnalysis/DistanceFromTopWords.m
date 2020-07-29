function DistanceFromTopWords(targetIdx,topNwords,topNcnts,libCodes,EcCount,varargin)
% 
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'N', 'positive', 1};
defaults(end+1,:) = {'targetIdx','positive',113:140};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A MList is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

k = 0;
blanks = zeros(28,1);
ctalk = zeros(28,1);
for j=targetIdx
h = zeros(parameters.N,1); 
for i=1:parameters.N
h(i) = sum(abs(topNwords(i,:) - libCodes(j,:)));
end
k=k+1;
blanks(k) = min(h);
ctalk(k) = sum(topNcnts.^(1./h));
end

figure(15); plot(blanks,EcCount(parameters.targetIdx),'k.'); 
xlabel('min hamming distance to top 15 words');
ylabel('counts'); 
title('L4E1 blank and Ecoli bits');
PresentationPlot('MarkerWidth',20);

figure(16); plot(ctalk,EcCount(parameters.targetIdx),'k.'); 
xlabel('"expected" cross-talk   sum(topNcnts.^(1./h)) ','Interpreter','none');
ylabel('counts'); 
title('L4E1 blank and Ecoli bits');
PresentationPlot('MarkerWidth',20);