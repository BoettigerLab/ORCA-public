function [finalHammingCodewords,finalSECDEDcodewords] = GetHammingCodebook()
% 
% 
% M = 4 allows 700+ SECDED words and 400 Hamming words, in 16 or 15
% readouts respectively, with max numOnBits (8 or 7 respectively)
% M = 4 with a fixed numOnBits of 4 is the SECDED-16-bit, 140 word code
% 
% M=5 is the 32-bit SECDED

% 
% N=2^M-1 where M=N-K.
% Balanced-ON Hamming codes at M=3 are too small
% pre-compute the M=3,4,5
M=4
N = 2^M-1
K = N-M
uncodedwords = rem(floor([0:2^K-1]'*pow2(-(K-1):0)),2);
hammingCodewords = encode(uncodedwords,N,K,'hamming');

numOnBits = []; % [];
figure(1); clf; hist(sum(hammingCodewords,2),1:N);
wordsPerOnBitCnt = hist(sum(hammingCodewords,2),1:N);
if isempty(numOnBits)
   [totMsgs, numOnBits] = max(wordsPerOnBitCnt);
end

finalHammingCodewords = hammingCodewords(sum(hammingCodewords,2)==numOnBits,:);
figure(1); clf; imagesc(finalHammingCodewords);  


%% SECDED
numOnBits = [];
figure(1); clf; hist(sum(hammingCodewords,2),1:N);
wordsPerOnBitCnt = hist(sum(hammingCodewords,2),1:N);
if isempty(numOnBits)
   [totMsgs, numOnBits] = max(wordsPerOnBitCnt);
end
inputWords = sum(hammingCodewords,2)==numOnBits |  sum(hammingCodewords,2)==numOnBits-1;
inputSecdedWords = hammingCodewords(inputWords,:);
numIn = size(inputSecdedWords,1);

SECDEDcodewords = zeros(numIn,N+1);
for n = 1:numIn
    SECDEDcodewords(n,:) = [mod(sum(hammingCodewords(n,:), 2), 2) hammingCodewords(n,:)]; %Compute SECDED codewords
end
figure(1); clf; hist(sum(SECDEDcodewords,2),1:N+1);
finalSECDEDcodewords = SECDEDcodewords(sum(SECDEDcodewords,2)==numOnBits+1,:);
figure(1); clf; imagesc(SECDEDcodewords); 

% figure(1); clf; imagesc(finalSECDEDcodewords); 