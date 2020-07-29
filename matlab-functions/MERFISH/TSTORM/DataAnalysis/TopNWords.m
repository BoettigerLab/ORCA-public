function [topNwords,topNcnts] = TopNWords(Mcents,libGenes,libCodes,varargin)
% Currently hardcoded N = 20 

N = 20;

[v,n] = occurrences(Mcents) ;
[nSort,idx] = sort(n);
topNwords = flipud(v(idx(end-N+1:end),:));
topNcnts = flipud(nSort(end-N+1:end)');

disp(['top ',num2str(N),' words:']);
disp(topNwords);

disp(['Counts for top ',num2str(N),' words:']);
disp(topNcnts); 

bookWords = ismember(libCodes,topNwords,'rows');
bookWordRanks = find(ismember(topNwords,libCodes,'rows'));
disp(['Codebook members in the top ',num2str(N),' words']);
disp(libGenes(bookWords))
disp(['rank of codebook members in top ',num2str(N),' words:',]);
disp(bookWordRanks)

