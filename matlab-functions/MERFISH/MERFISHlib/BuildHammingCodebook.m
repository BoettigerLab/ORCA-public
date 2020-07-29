function usedSECDEDcodewords = BuildHammingCodebook(numWords,varargin)
%
% 

% nMsg = 2^numDataBits; %Number of messages from the number of data bits.
if numWords < 14 + 1
    numBits = 8;
    numDataBits = 4;
elseif numWords < 38  +1  % approx'd
    numBits = 12;
    numDataBits = 11; % approx'd
elseif numWords < 140 + 1
    numBits = 16; %Total number of bits for SECDED codewords   
    numDataBits = 11; %Number of data bits
elseif numWords < 1001
    
end

defaults = cell(0,3);
defaults(end+1,:) = {'numBits','positive',[]}; % total number of bits (numColors x numHybes) that you intend to perform for this experiment
defaults(end+1,:) = {'onBits','positive',4}; % number of times each RNA will be labeled/detected. Equivelantly: the number of "1s" in each codeword 
defaults(end+1,:) = {'verbose','boolean',true}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);
   

onBits = pars.onBits; 


%% Main Function
%-------------------------------------------------------------------------

if numBits ~= 12
    % numDataBits = 11; %Number of data bits.
    uncodedwords = rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2); %Generate a matrix of uncoded words where the rows are the uncoded words
    hammingcodewords = encode(uncodedwords,numBits-1,numDataBits,'hamming/binary'); %Compute Hamming codewords
    SECDEDcodewords = zeros(size(hammingcodewords,1),numBits);
    for n = 1:length(hammingcodewords)
        SECDEDcodewords(n,:) = [mod(sum(hammingcodewords(n,:), 2), 2) hammingcodewords(n,:)]; %Compute SECDED codewords
    end               
    usedSECDEDcodewords = SECDEDcodewords(sum(SECDEDcodewords,2)==onBits,:);

%% 

else
% Make 12-Bit SECDED with 4 on bits (38 words)
%    hacked from a 16 bit SECDED

    numDataBits = 11;
    numBits = 16;
    uncodedwords = rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2); %Generate a matrix of uncoded words where the rows are the uncoded words
    hammingcodewords = encode(uncodedwords,numBits-1,numDataBits,'hamming/binary'); %Compute Hamming codewords
    SECDEDcodewords = zeros(size(hammingcodewords,1),numBits);
    for n = 1:length(hammingcodewords)
        SECDEDcodewords(n,:) = [mod(sum(hammingcodewords(n,:), 2), 2) hammingcodewords(n,:)]; %Compute SECDED codewords
    end 

    twelveBitSECDED = SECDEDcodewords(1:16:end,1:12);
    [nWords, nLetters] = size(twelveBitSECDED);
    hD_short = zeros(nWords);
    for j=1:nWords
    for i=1:nWords  % i = 2;
        hD_short(i,j) = sum(xor(twelveBitSECDED(i,:), twelveBitSECDED(j,:)));
    end
    end

    hD_short(hD_short==0) = inf;
    min(hD_short);
    bits12on4 = twelveBitSECDED(sum(twelveBitSECDED,2)==onBits,:);

    [nWords, nLetters] = size(bits12on4);
    hD_4onbits = zeros(nWords); 
    for j=1:nWords
    for i=1:nWords  % i = 2;
        hD_4onbits(i,j) = sum(xor(bits12on4(i,:), bits12on4(j,:)));
    end
    end
    % figure(2); clf; imagesc(hD_4onbits); colorbar;
    usedSECDEDcodewords = bits12on4;
end


