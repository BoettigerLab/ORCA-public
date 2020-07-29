function [codebook, code2ID] = Generate2NCodebook(N)

initialWords = dec2bin(0:(2^N-1));
codebook = cellstr(initialWords);
code2ID = containers.Map(codebook, 1:length(codebook));

