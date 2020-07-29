function [allCnts,geneNames] = ExtractCounts(codebook)

allCnts = cellfun(@Gcnts, {codebook.Sequence});
geneNames = cellfun(@Gname, {codebook.Sequence},'UniformOutput',false);

function cnts = Gcnts(x)
gaps = strfind(x,' ');
cnts = str2double(x(gaps(2):end));

function geneNames = Gname(x)
gaps = strfind(x,' ');
geneNames =  x(1:gaps(1)-1) ;