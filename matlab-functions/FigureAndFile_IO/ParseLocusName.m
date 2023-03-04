function [chr,locusStart,locusEnd] = ParseLocusName(locustxt,varargin)
%--------------------------------------------------------------------------
% parse a UCSC formatted chromatin locus into a chromosome name (string), 
% start coordinate (double), and end coordinate (double)
% 
%--------------------------------------------------------------------------
% Example
% [chr,locusStart,locusEnd] = ParseLocusName(chr2R:12893885-12993885)
% chr = chr2R, 
% locusStart = 12893885,
% locusEnd = 12993885,
% 
% See also WriteLocusName
%--------------------------------------------------------------------------

defaults = cell(0,3);
defaults(end+1,:) = {'removeChr','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

cellPassed = false; 

if iscell(locustxt)
    numLoci = length(locustxt);
    cellPassed = true; 
else
    numLoci = 1;
    locustxt = {locustxt};
end

locusStart = zeros(numLoci,1);
locusEnd = zeros(numLoci,1);
chr = cell(numLoci,1); 

for i = 1:numLoci
    try
     locusname = regexprep(locustxt{i},' ',''); % remove spaces 
     s = regexp(locusname,':');
     if ~isempty(s)
         sep = regexp(locusname,'-|\.\.');
         if pars.removeChr
            chr{i} = regexprep(locusname(1:s-1),{'chr','Chr'},{'',''});
         else
             chr{i} = locusname(1:s-1);
         end
         locusStart(i) = str2double(locusname(s+1:sep-1));
         locusEnd(i) = str2double(locusname(sep+1:end));
     else
         nameparts = strsplit(locusname,',');
         chr{i} = nameparts{1};
         locusStart(i) = str2double(nameparts{2});
         locusEnd(i) = str2double( nameparts{3});
     end
    catch er
        warning(er.message);
        error(['unable to parse locus name ', locustxt{i}]);
    end
end

if ~cellPassed
    chr = chr{1};
end