function [PmCount,EcCount] = PlotCountsPerGene(codebook,hammingDis,freqUniqueMsg,varargin)

%% Defaults
libGenes = '';

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'libGenes'
                libGenes = CheckParameter(parameterValue,'cell','libGenes');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

if isempty(libGenes)
libGenes = {codebook.Sequence};
end

%% Main Function
numGenes = length(codebook);
EcCount = zeros(numGenes,1);
PmCount = zeros(numGenes,1);
for n=1:numGenes
    PmCount(n) = sum(freqUniqueMsg(hammingDis(:,n)==0));
    EcCount(n) = sum(freqUniqueMsg(hammingDis(:,n)==1)) + ...
                     sum(freqUniqueMsg(hammingDis(:,n)==0));
end

numMatch = num2str(sum(PmCount));
numDecoded = num2str(sum(EcCount));
disp(['Perfect Matches ',numMatch]);
disp(['Clusters Decoded ',numDecoded]);

subplot(2,1,1);
bar(EcCount,'w');
hold on
geneNames = cell(numGenes,1); 
x = 1:numGenes;
for i=1:numGenes   
    geneNames{i} = libGenes{i};
    sp = strfind(libGenes{i},' ');
    if ~isempty(sp)
        geneNames{i} = geneNames{i}(1:sp(1)); 
    end  
    uLabel = text(x(i),1,['{\color{red}',geneNames{i},'}', ' ',...
                          '{\color{black}',num2str(EcCount(i)),'}'],...
                          'FontWeight','bold','FontName','FixedWidth');
    set(uLabel,'rotation',90)
end
title(['Error Corrected mRNA Counts per gene (',numDecoded,' clusters)']);
ylim([0,max( [max(EcCount),40]) ]);

subplot(2,1,2);
bar(PmCount,'w');
hold on
x = 1:numGenes;
for i=1:numGenes  
    uLabel = text(x(i),1,['{\color{red}',geneNames{i},'}', ' ',...
                          '{\color{black}',regexprep(codebook(i).Header,' ',''),'}'],...
                          'FontWeight','bold','FontName','FixedWidth');
    set(uLabel,'rotation',90)
end
ylim([0,max( [max(EcCount),40]) ]);
title(['Perfect match mRNA Counts per gene (',numMatch,' matches)']);
set(gcf,'color','w');

