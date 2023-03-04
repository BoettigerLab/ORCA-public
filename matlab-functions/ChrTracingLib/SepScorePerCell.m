function sepScores = SepScorePerCell(mapIn,b1,b2,varargin)
% single cell inter-domain vs. intra domain median ratio, for a given stack
% of distance maps "mapIn" in two specified domains b1 and b2.
%   intra-domain distance is the median of the b1 and b2 intradomain
%   inter-domain distance is the b1 vs b2.
%   score is inter/intra -- for distances, >1 is separate



defaults = cell(0,3);
defaults(end+1,:) = {'badHybes','integer',[]};
defaults(end+1,:) = {'normDistance','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[nB,~,nCells] = size(mapIn);
if pars.normDistance
    normMap = NormMap(nanmedian(mapIn,3));
else
    normMap = ones(nB);
end

sepScores = nan(nCells,1);
for c=1:nCells
    a = mapIn(:,:,c)./normMap;
    a(pars.badHybes,:) = nan;
    a(:,pars.badHybes) = nan;    
    L  = a(b1,b1);
    L = triu(L,1); 
    L(L==0) = nan;
    R = a(b2,b2);
    R = triu(R,1);
    R(R==0) = nan;
    X = a(b1,b2);
    X = tril(X,-1);
    X(X==0) = nan;
    inter = nanmedian(X(:));
    intra = nanmedian([L(:);R(:)]);
    sepScores(c) = inter./intra;
end