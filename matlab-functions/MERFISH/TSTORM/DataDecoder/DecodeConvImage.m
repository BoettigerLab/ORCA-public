function [perCnt,allCnt,spotMatrix,spotIDs,spotDaxMatrix,allCorCodes] = DecodeConvImage(wordsDetected,codebook,varargin)
% compare the codebook to the words read out at each mRNA centroid found in
% the image.  
% 4 ON-bit codebooks are hardcoded into this for speed and simplicity. 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showplots', 'boolean', false};
defaults(end+1,:) = {'brightnessPerSpot', 'array', []};
defaults(end+1,:) = {'daxValuePerSpot', 'array', []};
defaults(end+1,:) = {'correctErrors', 'boolean', true};
defaults(end+1,:) = {'libCodes', 'array', []};
defaults(end+1,:) = {'allCorCodes','cell',{}};
defaults(end+1,:) = {'codeType',{'binary','decimal'},'binary'};
defaults(end+1,:) = {'numHybes','positive',[]};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires flists, mRNAcents');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

showplots = parameters.showplots;
brightnessPerSpot = parameters.brightnessPerSpot; 
allCorCodes = parameters.allCorCodes;
daxValuePerSpot = parameters.daxValuePerSpot; 

    
% % libCodes may be prepassed for speed 
% if isempty(parameters.libCodes) && isstruct(codebook)
%     codebook = cellfun(@str2num,{codebook.Header},'UniformOutput',false);
%     codebook = cat(1,codebook{:});
% elseif isempty(parameters.libCodes) && strcmp(parameters.codeType,'binary'); 
%     codebook = codebook;
% elseif ~isempty(parameters.libCodes)
%     codebook = parameters.libCodes;
% end
% 


%% Main Function
if strcmp(parameters.codeType,'binary')
    decWords = bi2de(wordsDetected); % this type of indexing is fastest
    decCodes = bi2de(codebook);
    [numGenes,numHybes] = size(codebook);
elseif strcmp(parameters.codeType,'decimal') && parameters.correctErrors && isempty(allCorCodes)
    decWords = wordsDetected; % this type of indexing is fastest
    decCodes = codebook;
    if isempty(parameters.numHybes)
        error('either numHybes or allCorCodes is required for error correcting decimal codebooks'); 
    end
    codebook = de2bi(codebook,parameters.numHybes);
    [numGenes,numHybes] = size(codebook);
elseif strcmp(parameters.codeType,'decimal')
    decWords = wordsDetected; % this type of indexing is fastest
    decCodes = codebook;
    numGenes = length(decCodes); 
end
    

if isempty(brightnessPerSpot) && nargout > 2
    brightnessPerSpot = 1E3*ones(length(wordsDetected),numHybes);
end


% need
% 1). codebook BINARY 


% compute all word variants (gain or loss of 1 bit)
if parameters.correctErrors && isempty(allCorCodes); 
    blankCode = false(1,numHybes); 
    allCorCodes = cell(numGenes,1); % all correctable codes
    for n=1:numGenes
        % codes(n,:) = str2num(codebook(n).Header);
        allCodesMatrix = zeros(numHybes*2,numHybes);
        for h=1:numHybes
            gainBit = blankCode; 
            gainBit(h) = true; 
            lossBit = ~blankCode;
            lossBit(h) = false; 
            allCodesMatrix(2*(h-1)+1:2*(h-1)+2,:) = [codebook(n,:) | gainBit; codebook(n,:) & lossBit ];
        end
        correctableCodes = sum(allCodesMatrix,2)==3 | sum(allCodesMatrix,2)==5;
        allCorCodes{n} = bi2de(allCodesMatrix( correctableCodes, :));
    end
end

% Decode all words in wordsDetected Matrix
%   accomplished by only looping over genes, not total words detected.  
spotMatrix = cell(numGenes,1); 
spotDaxMatrix = cell(numGenes,1); 
spotPerIDs = zeros(length(wordsDetected),1);
spotIDs  = zeros(length(wordsDetected),1);
perCnt = zeros(numGenes,1); 
corCnt = zeros(numGenes,1);
for n = 1:numGenes   % n =9
    perMatches = decWords==decCodes(n);   
    perCnt(n) = sum(perMatches);
    
    if nargout > 2
        spotPerIDs(perMatches) = n; 
        spotIDs(perMatches) = n; 
    end
    
    if parameters.correctErrors
        [~,allErrMatches] = intersect(decWords,allCorCodes{n});
        corCnt(n) = length(allErrMatches);
        spotIDs(allErrMatches) = n; 
    end
    
    if ~isempty(brightnessPerSpot)
        spotMatrix{n} = brightnessPerSpot(spotIDs==n,:)';
    end
    if ~isempty(daxValuePerSpot)
        spotDaxMatrix{n}= daxValuePerSpot(spotIDs==n,:)';
    end
end 

allCnt = perCnt+corCnt;

if ~parameters.correctErrors
    spotMatrix = wordsDetected;
end

% n=3;
% figure(1); clf; imagesc(  spotMatrix{n} ); 
% colormap(jet(256));
% xlabel('mRNA #'); ylabel('hybe');
% title(libGenes{n});
% PresentationPlot();

