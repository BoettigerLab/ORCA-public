function codebook = BuildNPCode(numDataBits,varargin)
%  for a given number of barcodes/databits N, return the largest NchooseP
%  codebook
% P (numOnBits) can either be selected to maximize the codebook length or
% be a fixed parameter passed by the user.
% 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'numOnBits','integer',0}; % 0 will return num OnBits that gives the largest codebook
defaults(end+1,:) = {'showOnBits','integer',1};
defaults(end+1,:) = {'showCodebook','integer',2};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% numDataBits = 9;

%Generate a matrix of uncoded words where the rows are the uncoded words
uncodedwords = rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2);

% determine the number of On Bits that produces the max number of codes
numOnBits = pars.numOnBits; % % []; % [];
if pars.showOnBits
figure(1); clf; hist(sum(uncodedwords,2),1:N);
end
wordsPerOnBitCnt = hist(sum(uncodedwords,2),1:N);
if numOnBits == 0
[totMsgs, numOnBits] = max(wordsPerOnBitCnt);
end

% keep only combinations with the requested (or max) number of on bits
codebook = uncodedwords(sum(uncodedwords,2)==numOnBits,:);
if pars.showCodebook
figure(pars.showCodebook); clf; imagesc(codebook);  
end
   