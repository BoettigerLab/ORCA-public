function [matchedbases,deletedbases,insertedbases,otherbases] = ...
    parseCIGAR(cigar) 
%--------------------------------------------------------------------------
% Parse a CIGAR string from a SAM file 
% 
%--------------------------------------------------------------------------% 
% Inputs:
% 
% cigar / string -- a CIGAR string from a SAM file.  can be retrieved from
%                   a bioIndexed SAM file by readData = read(bioIndexed)
%                   cigar = readData.CigarString.
%
%--------------------------------------------------------------------------
% Outputs
% 
% matchedbases / double  -- total number of matched or mismatched bases 
% deletedbases / double  -- total deleted bases in string
% insertedbases / double -- total number of inserted bases
% otherbases / double -- total number of other annotations (usually S). 
%
%--------------------------------------------------------------------------
% Example 1:
% cigar = '5S31M1I57M2D20S'
% matchedbases = 31 + 57
% deletedbases = 2
% insertedbases = 1
% otherbases = 25
% 
% Example 2:
% probeCIGARs is a cell array of CIGAR strings
% [baseMatches,baseDeletions] = cellfun(@parseCIGAR,probeCIGARs,...
%                                     'UniformOutput',false);
% baseMatches = cell2mat(baseMatches);
% baseDeletions = cell2mat(baseDeletions);
% vectors with the number of base matches and base deletions recorded in
% each CIGAR string.  
%--------------------------------------------------------------------------


%% Main Function

cidx = [0,regexp(cigar,'[A-Z]')];
matchedbases = 0;
deletedbases = 0;
insertedbases = 0;
otherbases = 0; 
for c=1:length(cidx)-1 % c = 2
    matchtype = cigar(cidx(c+1));
    matchbases = str2double(cigar(1+cidx(c):cidx(c+1)-1));
    if strcmp(matchtype,'M')
        matchedbases = matchedbases + matchbases;
    elseif strcmp(matchtype,'D')
        deletedbases = deletedbases + matchbases;
    elseif strcmp(matchtype,'I');
        insertedbases = insertedbases + matchbases;
    else
        otherbases = otherbases+matchbases;
    end
end
