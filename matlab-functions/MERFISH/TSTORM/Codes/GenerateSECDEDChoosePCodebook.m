function [exactCodebook, allCode2GeneID, correctedCodebook, exactCode2GeneID, correctedCode2GeneID] = GenerateSECDEDChoosePCodebook(numDataBits, P)

% Generate all words then strip to number of on bits
words = GenerateExtendedHammingWords(numDataBits);
words = words(sum(words,2)==P,:);

% Covert words to exact codebook
exactCodebook = ConvertLogicalWordsToString(words);

% Calculate correctedCodebook
correctedCodebook = {};
for i=1:size(words,1)
    correctedCodebook{i} = ConvertLogicalWordsToString(cell2mat(GenerateSurroundingCodewords(words(i,:)==1,1)));
end

% Create map function
exactCode2GeneID = containers.Map(cellstr(exactCodebook), 1:size(exactCodebook,1));

% Create map function
geneInds = repmat([1:size(exactCodebook,1)]', [1 size(words,2)])';
geneInds = reshape(geneInds, [numel(geneInds) 1]);
correctedCode2GeneID = containers.Map(cellstr(cat(1,correctedCodebook{:})), geneInds);

% Convert codebook to cell
exactCodebook = cellstr(exactCodebook);

% Create exact and corrected map
allCode2GeneID = [exactCode2GeneID; correctedCode2GeneID];

end


function stringWords = ConvertLogicalWordsToString(logicalWords)
    stringWords = repmat('', size(logicalWords));
    stringWords(find(logicalWords(:))) = '1';
    stringWords(find(~logicalWords(:))) = '0';
    stringWords = reshape(stringWords, size(logicalWords));
end