function wordsWithCorrBits = WordsWithCorrBits(bitPairs,libCodes,libGenes)

% bitPairs = [1,2; 1,3; 2,3; 2,4];
numPairs = size(bitPairs,1);

wordsWithCorrBits = cell(numPairs,2); 
for n=1:numPairs
idx = libCodes(:,bitPairs(n,1)) + libCodes(:,bitPairs(n,2)) == 2;
disp('Genes on during hybes '); disp(bitPairs(n,:));
disp(libGenes(idx));
wordsWithCorrBits{n} = {bitPairs(n,:), libGenes(idx)};
end