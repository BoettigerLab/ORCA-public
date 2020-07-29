function [probeData,libExpect,libExpectPerHybe,libGenes,libCodes] = ExpectedLocsPerHybe2(codebook,probeData,varargin)
% [probeData,libExpect,libExpectPerHybe,libGenes,libCodes] = ExpectedLocsPerHybe2(codebook,probeData,varargin)

% % Defaults
% fpkmPath = '\\cajal\TSTORMdata\GenomeData\A549\';
% codebookPath = '';
% 
% %--------------------------------------------------------------------------
% %% Parse variable input
% %--------------------------------------------------------------------------
% if nargin > 1
%     if (mod(length(varargin), 2) ~= 0 ),
%         error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
%     end
%     parameterCount = length(varargin)/2;
%     for parameterIndex = 1:parameterCount,
%         parameterName = varargin{parameterIndex*2 - 1};
%         parameterValue = varargin{parameterIndex*2};
%         switch parameterName
%             case 'fpkmPath'
%                 fpkmPath = CheckParameter(parameterValue,'string','fpkmPath');
%             case 'codebookPath'
%                 codebookPath = CheckParameter(parameterValue,'string','codebookPath');
%             otherwise
%                 error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
%         end
%     end
% end

%%


% Extract gene names from codebook
numHybes = length(regexprep(codebook(1).Header,' ','')); % number of stains 
libNames = {codebook.Sequence}';
gaps = cellfun(@(x) strfind(x,' '),libNames,'UniformOutput',false);
libGenes = cellfun(@(x,y) x(1:y(1)-1),libNames,gaps,'UniformOutput',false);
if length(gaps{1}) > 15
    libIsoIds = cellfun(@(x,y) x(y(4)+1:y(5)-1),libNames,gaps,'UniformOutput',false);
else
    libIsoIds = '';
end

libExpect = probeData.FPKM;
% setdiff(probeData.CommonName,libGenes)

libCodes = logical(cell2mat(cellfun(@str2num, {codebook.Header},'UniformOutput',false)'));
libExpectPerHybe = zeros(numHybes,1);
for i=1:numHybes % 
    libExpectPerHybe(i) = sum(libExpect(libCodes(:,i)));
end
libExpectPerHybe = libExpectPerHybe/sum(libExpectPerHybe);

probeData.codewords = cellstr(num2str(libCodes));

