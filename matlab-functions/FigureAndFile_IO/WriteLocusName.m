function locusname = WriteLocusName(chr,locusStart,locusEnd)
% returns formatted string locus name from numeric locus data. 
% if chr is just a number, assumes this is indexing the Drosophila
% chromosomes (2L, 2R, 3L, 3R, 4, X)

%% 

numNames = length(locusStart); 

%% single entry
if numNames == 1
    if ~ischar(chr)
        [~, ~, chrNames] = GetColorIntervals('showPlots',false);
        chr = chrNames{chr};
    end   
    locusname = [chr,':',num2str(round(locusStart)),'-',num2str(round(locusEnd))];

%% multi entry
elseif numNames > 1
    locusname = cell(numNames,1);
    for i=1:numNames
        if length(chr) > 1 && iscell(chr)
            chr_i = chr{i};
        elseif length(chr) > 1 && ~iscell(chr)
            chr_i = chr(i);
        else 
            chr_i = chr;
        end

        if ~ischar(chr_i )
            [~, ~, chrNames] = GetColorIntervals('showPlots',false);
            chr_i  = chrNames{chr_i };
        end

        locusname{i} = [chr_i,':',num2str(round(locusStart(i))),'-',num2str(round(locusEnd(i)))];
    end
end
 