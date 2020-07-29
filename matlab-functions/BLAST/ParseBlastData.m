function blastResults = ParseBlastData(blastData)
%--------------------------------------------------------------------------
%  blastResults = ParseBlastData(blastData)
% 
% parses the blast data structure produced by matlab into a structure of
% identical length indexed arrays, which is faster and much easier to
% navigate. 

%% 
numQueries = length(blastData);

blastResults.topHitName = cell(numQueries,1);
blastResults.secHitName = cell(numQueries,1);
blastResults.topHitLength = zeros(numQueries,1);
blastResults.secHitLength = zeros(numQueries,1);
blastResults.topHitBpsMatched = zeros(numQueries,1);
blastResults.secHitBpsMatched = zeros(numQueries,1);
blastResults.topQuerySeq = cell(numQueries,1);
blastResults.topHitSeq = cell(numQueries,1);
blastResults.secQuerySeq = cell(numQueries,1);
blastResults.secHitSeq = cell(numQueries,1);
blastResults.query = cell(numQueries,1);
blastResults.numHits = zeros(numQueries,1);  % 

for n=1:numQueries
    blastResults.query{n} = blastData(n).Query;
    blastResults.numHits(n) = length(blastData(n).Hits);
    hitTemp =  blastData(n).Hits;
    if ~isempty(hitTemp)
        blastResults.topHitName{n} = hitTemp(1).Name;
        blastResults.topHitLength(n) = hitTemp(1).Length;
        blastResults.topHitBpsMatched(n) = length(strfind(hitTemp(1).HSPs(1).Alignment(2,:),'|'));
        blastResults.topQuerySeq{n} = hitTemp(1).HSPs(1).Alignment(1,:);
        blastResults.topHitSeq{n} = hitTemp(1).HSPs(1).Alignment(3,:);
    end
    if length(hitTemp) > 1
        blastResults.secHitName{n} = hitTemp(2).Name;
        blastResults.secHitLength(n) = hitTemp(2).Length;
        blastResults.secHitBpsMatched(n) = length(strfind(hitTemp(1).HSPs(1).Alignment(2,:),'|'));
        blastResults.secQuerySeq{n} = hitTemp(2).HSPs(1).Alignment(1,:);
        blastResults.secHitSeq{n} = hitTemp(2).HSPs(1).Alignment(3,:);
    end
end
 
