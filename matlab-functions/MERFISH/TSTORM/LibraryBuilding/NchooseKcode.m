function codebook = NchooseKcode(numHybes,numOnBits)
% generates an n choose k codebook 

% numHybes = 14;
% numOnBits = 4;

nCk = nchoosek(1:numHybes,numOnBits);
numWords = size(nCk,1);
codebook = false(numWords,numHybes);
for i=1:numWords
    codebook(i,nCk(i,:)) = true;
end
