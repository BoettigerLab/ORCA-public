function parameters = BarGraphAllGenes(rnaCnts,libGenes,varargin) 
 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'logscale', 'boolean', false};
defaults(end+1,:) = {'codebook','struct',[]};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'counts and names are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%%

codebook = parameters.codebook;

 EcCount = sum(rnaCnts,2);
 numGenes = length(EcCount);
 
 negIdx = [StringFind(libGenes,'blank');StringFind(libGenes,'notarget')];
 negCnt = max(EcCount(negIdx));

if parameters.logscale 
    bar(log10(EcCount),'w');
else
    bar(EcCount,'w');
end

numDecoded = sum(EcCount);
hold on
geneNames = cell(numGenes,1); 
x = 1:numGenes;
for i=1:numGenes   
    geneNames{i} = libGenes{i};
    sp = strfind(libGenes{i},' ');
    if ~isempty(sp)
        geneNames{i} = geneNames{i}(1:sp(1)); 
    end  
    geneCnt = EcCount(i);
    if geneCnt > negCnt
        nameclr = 'blue';
    else
        nameclr = 'red';
    end
    if isempty(codebook)
        geneCode = '';
    else
        geneCode = ['   ','{\color{black}',regexprep(codebook(i).Header,' ',''),'}'];
    end
    
    uLabel = text(x(i),0,['{\color{',nameclr,'}',geneNames{i},'}', ' ',...
                          '{\color{red}',num2str(EcCount(i)),'}',...
                           geneCode],...
                          'FontWeight','bold','FontName','FixedWidth');
    set(uLabel,'rotation',90)
end
title(['Error Corrected mRNA Counts per gene (',num2str(numDecoded),' clusters)']);
if parameters.logscale
    ylim([0,max( [max(log10(EcCount)),0]) ]);
else
    ylim([0,max( [max(EcCount),0]) ]);
end
set(gcf,'color','w');