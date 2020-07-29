function [codebook, code2ID] = GenerateNChoosePCodebook(N, P)

initialWords = dec2bin(0:(2^N-1));
numOn = sum(initialWords == '1',2);

goodInds = numOn == P;

codebook = cellstr(initialWords(goodInds,:));

code2ID = containers.Map(codebook, 1:length(codebook));

