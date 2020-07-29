function corruptedWords = SimulateMeasurementError(actualWords, probability)

strings = '10';
corruptedWords = actualWords;

for i=1:size(actualWords,1)
    onInds = actualWords(i,:) == '1';
    localProb(onInds) = probability(1,onInds);
    localProb(~onInds) = probability(2,~onInds);
    
    bitsToFlip = rand(1, length(localProb)) <= localProb;
    
    previousValues = corruptedWords(i,bitsToFlip) == '1';
    newValues = previousValues + 1;
    newValues = strings(newValues);
    
    corruptedWords(i,bitsToFlip) = newValues;
end


