function data = BLASTreadLocal(blastreport, mopt)
%BLASTREADLOCAL read data from local BLAST report
%
%   DATA = BLASTREADLOCAL(BLASTREPORT, FORMAT) reads BLASTREPORT, a locally
%   created BLAST report text file, and returns DATA, a MATLAB structure or
%   array of structures (if multiple query sequences) containing fields
%   corresponding to BLAST keywords. BLASTREPORT can be a text file name or
%   a one-dimensional char array containing the text for a local BLAST
%   report. The number and type of fields in DATA depends on FORMAT, an
%   integer specifying the alignment format used to create BLASTREPORT.
%   Valid FORMAT options are:
%
%        0  - pairwise
%        1  - query-anchored, showing identities
%        2  - query-anchored, no identities
%        3  - flat query-anchored showing identities
%        4  - flat query-anchored, no identities
%        5  - query-anchored, no identities and blunt ends
%        6  - flat query-anchored, no identities and blunt ends
%        8  - tabular
%        9  - tabular with comment lines
%
%   Examples:
%
%   % Example 1
%   % Retrieve two sequences from GenBank, and write them into a fasta file,
%   % using the accession numbers as headers.
%   S1 = getgenbank('M28570.1');
%   S2 = getgenbank('M12565');
%   Seqs(1).Header = S1.Accession;
%   Seqs(1).Sequence = S1.Sequence;
%   Seqs(2).Header = S2.Accession;
%   Seqs(2).Sequence = S2.Sequence;
%   fastawrite('query_multi_nt.fa', Seqs);
% 
%   % Please see the examples in BLASTFORMAT help on how to create a
%   % local blastable database ecoli.nt.
% 
%   % Search query sequences against a local database (ecoli.nt), and save
%   % the results in a tabular report file.
%   blastlocal('inputquery', 'query_multi_nt.fa', ... 
%              'database', 'ecoli.nt', 'tofile', 'myecoli_nt8.txt', ...
%              'program', 'blastn', 'format', 8)
%
%   % Read BLAST report
%   blastreadlocal('myecoli_nt8.txt', 8)
%
%   % Example 2
%   % Search query sequences and write a query-anchored report file
%   blastlocal('inputquery', 'query_multi_nt.fa', ...
%              'database', 'ecoli.nt', 'tofile', 'myecoli_nt1.txt', ...
%              'program', 'blastn', 'format', 1);
%
%   % Read BLAST report
%   results = blastreadlocal('myecoli_nt1.txt', 1)
%
% See also BLASTFORMAT, BLASTLOCAL, BLASTNCBI, BLASTREAD, GETBLAST.

% References
% http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs

% Copyright 2007-2016 The MathWorks, Inc.


if nargin < 2
	error(message('bioinfo:blastreadlocal:NotEnoughInput'));
end

data = [];

if ~isnumeric(mopt)|| ~isscalar(mopt) || mopt == 7 || mopt > 9
	error(message('bioinfo:blastreadlocal:InvalidOption'));
end

if isempty(blastreport)
	warning(message('bioinfo:blastreadlocal:NoHits'));
	return
end

if ~ischar(blastreport)
	error(message('bioinfo:blastreadlocal:NotChar'));
else
	if isvector(blastreport)
		if size(blastreport,2) == 1
			blastreport = blastreport';
		end
	end
	
	if  exist(blastreport, 'file') % is it a file?
		try
			fid = fopen(blastreport);
			report = fread(fid,'*char')';
			fclose(fid);
		catch theErr
			if strcmpi(theErr.identifier,'MATLAB:nomem')
				error(message('bioinfo:blastreadlocal:ReportTooBigToOpen'));
			else
				error(message('bioinfo:blastreadlocal:CannotReadInput', blastreport));
			end
		end
	else % is it a file name or actual report?
		if isempty(regexpi(blastreport, '\n', 'once'))  % must be a file name
			error(message('bioinfo:blastreadlocal:CannotFindInput', blastreport));
		else % must be a one-dim char array with the actual report
			report = blastreport;
		end
	end
end

%=== main
try
    if mopt == 0
        data = blastParser0(report);
    elseif mopt == 8
        data = blastParser8(report);
    elseif mopt == 9
        alg = regexprep(regexpi(report, '#\s+((T*)BLAST[N|P|X]\s+\d+\.\d+\.\d+\s+\[\w{3}-\d{2}-\d{4}\])', 'match'), '^#\s+', '');
        databases = regexprep(regexpi(report, '# Database: (.*?)\s+', 'match'), '^# Database: ', '');
        databases = regexprep(databases, '[\r\n]+','');
        report =  regexprep(report, '#(.*?)bit score',''); % remove all comment lines
        data = blastParser8(report);
        if ~isempty(data)
            for i = 1:numel(alg)
                data(i).Algorithm = alg{i};
                data(i).Database = databases{i};
            end
            data = orderfields(data, {'Algorithm', 'Query', 'Database', 'Hits'});
        end
    else % mopt== 1 to 6
        data = blastParser1to6(report);
    end
catch ME
	if strcmpi(ME.identifier,'MATLAB:nomem')
		error(message('bioinfo:blastreadlocal:ReportTooBigForParsing'));
	elseif strncmpi(ME.identifier,'bioinfo:',8) % did we catch one of our errors?
		rethrow(ME)
	else
		error(message('bioinfo:blastreadlocal:UnableToRead'))
	end
end

if isempty(data)
	warning(message('bioinfo:blastreadlocal:NoHits'));
end


%==========================================================================
% SUBFUNCTIONS
%==========================================================================
function blasttext = removeHTMLTags(blasttext)
% Removes all HTML tags

% Find where the PREformatted report starts, usually the <PRE> tag followed
% by a BLAST program line:

s = regexp(blasttext,'<PRE>\s+(T*)BLAST[N|P|X]\s+\d+\.\d+\.\d+\s{0,3}([\[\(]\w{3}-\d{2}-\d{4}[\]\)])?');

if isempty(s)
    % If not found, just try to find the <PRE> tag:
    s =  regexp(blasttext,'<PRE>');
end

if numel(s)>1
   warning(message('bioinfo:blastreadlocal:MultiplePreformattedSections'));
   s = s(1);
end

if ~isempty(s)
    blasttext = blasttext(s+5:end); % start after the <PRE> tag
end

% Removes all open/close HTML tags and the text between them
blasttext = regexprep(blasttext,'<(\w+).*?>.*?</\1>','');

% Remove remaining single HTML tags
blasttext = regexprep(blasttext,'<.*>','');

%==========================================================================
function [a, blasttext] = processHeader(blasttext)
% Removes all header information for each query. A header is defined as all
% the lines between the BLAST line and the start of the query (query=). It
% also figures out the BLAST program used.

[t, a] = regexpi(blasttext, ...
	'((T*)BLAST[N|P|X]\s+\d+\.\d+\.\d+\s{0,3}([\[\(]\w{3}-\d{2}-\d{4}[\]\)])?)', 'start', 'match'); % title starting line and algorithm

% A report can have the BLAST line several times, if so we check that they
% belong to the same program, version and date (if present), if not we
% warn. We keep this info once.
if numel(a)>1
    if numel(unique(a))>1
        warning(message('bioinfo:blastreadlocal:MultipleAlgorithms'));
    end
    a = a(1);
end

if isempty(a)
    warning(message('bioinfo:blastreadlocal:UnknownAlgorithm'));
    a = {'Unknown program'};
    return % nothing needs to be trimmed
end

% Find where the qeuries start
q = regexpi(blasttext,'Query= ');

% For each BLAST line fine its closest query start
z = true(size(blasttext));
for i = 1:numel(t)
    x = q(find(q>t(i),1))-1;
    z(t(i):x) = false;
end

% Trim headers
blasttext = blasttext(z);

%==========================================================================
function [d, blasttext] = processDatabase(blasttext)
% Removes all database information from blasttext. Database information can
% be stored anywhere: within the queries, in the header, or within the end
% stats. Database info has one of the following two formats:
%
% 1. 
%
%  Database: /Users/any/mydatabase.fa
%    Posted date:  Jun 1, 2010  11:00 PM
%  Number of letters in database: 9,752,417
%  Number of sequences in database:  10
%
% 2. 
%
%Database: ecoli.nt 
%400 sequences; 4,662,239 total letters
%

% Search first format:
[s,e,t] = regexp(blasttext, ['[\r\n]  Database:\s+(.*?)' ...
                            '    Posted date:\s+(.*?)'...
                            '  Number of letters in database:\s+([\d,]+)\s+'...
                            '  Number of sequences in database:\s+([\d,]+)'], ...
                            'start','end','tokens');

if isempty(t)
    dbNames = {};
    dbDates = {};
    dbNuLet = {};
    dbNuSeq = {};    
else
    t = [t{:}];
    dbNames = t(1:4:end);
    dbDates = t(2:4:end);
    dbNuLet = t(3:4:end);
    dbNuSeq = t(4:4:end);
end

% Remove database sections that match the first format
z = true(size(blasttext));
for i = 1:numel(s)
    z(s(i)+1:e(i)) = false;
end
blasttext = blasttext(z);

% Search second format:
[s,e,t] = regexp(blasttext, ['[\r\n]Database:\s+(.*?)' ...
                             '\s+([\d,]+)\s+sequences;'...
                             '\s+([\d,]+)\s+total letters'],...
                             'start','end','tokens');
                         
if ~isempty(t)
    n = numel(t);
    t = [t{:}]; 
    dbNames = [dbNames t(1:3:end)];
    dbDates = [repmat({''},1,n) t(2:3:end)];
    dbNuLet = [dbNuLet t(3:3:end)];
    dbNuSeq = [dbNuSeq t(2:3:end)];
end

% Remove database sections that match the second format
z = true(size(blasttext));
for i = 1:numel(s)
    z(s(i)+1:e(i)) = false;
end
blasttext = blasttext(z);

if numel(dbNames)>1
    if numel(unique(regexprep(dbNames,'\s','')))>1
       warning(message('bioinfo:blastreadlocal:MultipleDataBases'));
    end
end

if isempty(dbNames)
    warning(message('bioinfo:blastreadlocal:UnknownDataBase'));
    d.Name = {'Unknown database'};
    d.Date = '';
    d.NumberOfLetters = [];
    d.NumberOfSequences = [];
    return % nothing needs to be trimmed
end

d.Name = strtrim(regexprep(dbNames{1},'\s+',' '));
d.Date = strtrim(dbDates{1});
d.NumberOfLetters = str2double(dbNuLet{1});
d.NumberOfSequences = str2double(dbNuSeq{1});

%==========================================================================
function [st,blasttext] = processStats(blasttext)
% Parses and removes all the information relevant to stats. Stats may
% appear only at the end of the file or for every query. Either case, they
% always start with the 'Lambda K H' string. If they appear multiple times,
% the last appearance may have also some additional global stats. Since
% BLASTLOCALREAD returns the Statistics field for each query, the global
% stats are replicated into the query local stats.

[s,e,st] = regexp(blasttext, ['Lambda\s+K\s+H'...
                             '\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+'...
                             '(\s+Gapped\s+Lambda\s+K\s+H'...
                              '\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+)?'...
                             '(\s+Effective search space used:\s+[\d]+)?'],...
                             'start','end','match');
if isempty(st)
    warning(message('bioinfo:blastreadlocal:UnknownStatistics'));
    st = {'Unknown statistics'};
    return % nothing needs to be trimmed
end
                      
st = strcat(st,blasttext(e(end)+1:end));
st = regexprep(st,'[\n\r]+','\n');

% Remove database sections that match the stats format
e(end) = numel(blasttext);
z = true(size(blasttext));
for i = 1:numel(s)
    z(s(i):e(i)) = false;
end
blasttext = blasttext(z);

%==========================================================================
function [queryNames,queryLengths,q] = processQueries(blasttext)
% Find query regions and parse name and length, two formats are allowed:
%
% 1.
% Query= query_name   Length= number
%
% 2.
% Query= query_name (number letters)

[q,t] = regexpi(blasttext,'Query= (.*?)Length=(\d+)','start','tokens');
if isempty(q)
    [q,t] = regexpi(blasttext,'Query= (.*?)\(\s*(\d+)\s+letters\s*\)','start','tokens');
end

if isempty(q)
    error(message('bioinfo:blastreadlocal:QueryNotFound'));
end

t = [t{:}]; 
queryNames = strtrim(regexprep(t(1:2:end),'\s+',' '));
queryLengths = str2double(t(2:2:end));
q = [q numel(blasttext)+1]; % update index to queries to use it inside a loop for all the queries

%==========================================================================
function summInfo = summaryParser(summ)
% Extract information from the summary section: subject names, subject
% descriptions, scores and e-values.

s = regexpi(summ,'(?<name>([\w|\.]+))\s+(?<desc>([\.-\w\s]*))\s+(?<score>(\d+))\s+(?<evalue>([e-\d\.]*))', 'names');
for j = 1:numel(s) % for each subject j in the summary
	%summInfo(j).Name = s(j).name;
	%summInfo(j).Description = s(j).desc;
	summInfo(j).Name = [s(j).name ' ' s(j).desc]; %#ok<AGROW>
	summInfo(j).Score = str2double(s(j).score); %#ok<AGROW>
	summInfo(j).Expect = str2double(regexprep(s(j).evalue, '^(e[-+])',' 1$1')); %#ok<AGROW> % correct for abbreviated evalues
end

%==========================================================================
function out = blastParser1to6(blasttext)
% Parse all types of query-anchored blast reports (produced using local
% BLAST options m = 1-6).

%=== Remove all HTML formatted text and tags (if any)
blasttext = removeHTMLTags(blasttext);

%=== Search for databases in the BLAST report, parse them and remove them 
%    from the text:
[dataBase, blasttext] = processDatabase(blasttext);

%=== Search for the header of every query in the BLAST report, parse the
%    algorithm name and and remove them from the text:
[alg, blasttext] = processHeader(blasttext);

%=== Search for regions with statistics (either for each query or global
%    statistics), parse the regions and remove them from the text:
[stats,blasttext] = processStats(blasttext);

%=== Now, we should only have clean "Query=" regions, with/without
%    summary, with zero, one or more hits and each hit with one or more
%    HSPs.

%=== Parse the query names and their lengths:
[queryNames,queryLengths,q] = processQueries(blasttext);
numQueries = numel(queryNames);

if (numQueries~=numel(stats)) && (numel(stats)>1)
    warning(message('bioinfo:blastreadlocal:MultipleStatistics'));
    stats = stats(1);  %Ignoring non-matching stats and using only the first set found.  
end

%=== Preallocate space for each query
dummy1.Algorithm = '';
dummy1.Query = '';
dummy1.Length = [];
dummy1.Database = '';
dummy1.Hits.Name = '';
dummy1.Hits.Score = [];
dummy1.Hits.Expect = [];
dummy1.Alignment = [];
dummy1.Statistics ='';
out = repmat(dummy1,1,numQueries);

%=== Parse each query report
for i = 1:numQueries
    out(i).Algorithm = alg{1};
    out(i).Query = queryNames{i};
    out(i).Length = queryLengths(i);
    out(i).Database = dataBase.Name;
    
    queryReport = blasttext(q(i):(q(i+1)-1));
    
    if ~contains(queryReport,'No hits found *****')
        
        summStarts = strfind(queryReport,'Sequences producing significant alignments:');
        summEnd = regexpi(queryReport, '\n(\d+_0)','once'); %query label used in the alignments
        
        %=== parse summary
        summ = queryReport(summStarts:(summEnd-1));
        out(i).Hits = summaryParser(summ);
        
        %=== copy over the alignment info
        out(i).Alignment = regexprep(queryReport(summEnd:end),'[\n\r]+','\n');
        
    end % No hits found
    
    if numel(stats)==1
        out(i).Statistics = stats{1};
    else
        out(i).Statistics = stats{i};
    end
end % for i       
        
%==========================================================================
function out = blastParser0(blasttext)
% Parse pairwise BLAST reports (produced using local BLAST option m = 0,
% default).

%=== Remove all HTML formatted text and tags (if any)
blasttext = removeHTMLTags(blasttext);

%=== Search for databases in the BLAST report, parse them and remove them 
%    from the text:
[dataBase, blasttext] = processDatabase(blasttext);

%=== Search for the header of every query in the BLAST report, parse the
%    algorithm name and and remove them from the text:
[alg, blasttext] = processHeader(blasttext);

%=== Search for regions with statistics (either for each query or global
%    statistics), parse the regions and remove them from the text:
[stats,blasttext] = processStats(blasttext);

%=== Now, we should only have clean "Query=" regions, with/without
%    summary, with zero, one or more hits and each hit with one or more
%    HSPs.

%=== Parse the query names and their lengths:
[queryNames,queryLengths,q] = processQueries(blasttext);
numQueries = numel(queryNames);

if (numQueries~=numel(stats)) && (numel(stats)>1)
    warning(message('bioinfo:blastreadlocal:MultipleStatistics'));
    stats = stats(1);  %Ignoring non-matching stats and using only the first set found.  
end

%=== Determine type of BLAST search
use_frame = false;
use_positives = false;
use_strand = false;
if contains(blasttext,'Frame =') % Frame is included in translated searches
    use_frame = true;
end
if contains(blasttext,'Positives =') % Positives are included with peptide sequences
    use_positives = true;
end
if contains(blasttext,'Strand =') % Strand is included with nucleotide sequences
    use_strand = true;
end

maxHSPs = 100;

%=== Preallocate space for each query
dummy1.Algorithm = '';
dummy1.Query = '';
dummy1.Length = [];
dummy1.Database = '';
dummy1.Hits= [];
dummy1.Statistics ='';
out = repmat(dummy1,1,numQueries);
times = zeros(numQueries,1);
%=== Parse each query report
for i = 1:numQueries  % i=36
    tic;
    out(i).Algorithm = alg{1};
    out(i).Query = queryNames{i};
    out(i).Length = queryLengths(i);
    out(i).Database = dataBase.Name;
    
    queryReport = blasttext(q(i):q(i+1)-1);
    
    if ~contains(queryReport,'No hits found *****')
        
        %=== Extract subjects (hits) info relative to query i
        subjectStarts = regexpi(queryReport, '\n>');
        numSubjects = length(subjectStarts);
        if numSubjects > 1
            subjectEnds = [subjectStarts(2:end)-1 numel(queryReport)];
        else
            subjectEnds = numel(queryReport);
        end
        %s = regexpi(queryReport,'>(?<subject>([\w.|]+))\s+(?<desc>([\w\s-]*))\s+Length\s*=\s*(?<length>(\d+))', 'names');
        s = regexpi(queryReport, '>(?<subject>(.*?))\s+Length\s*=\s*(?<length>(\d+))', 'names');
        
        %=== Preallocate space for each subject
        dummy2.Name = '';
        dummy2.Length = [];
        dummy2.HSPs = [];
        r = repmat(dummy2, 1, numSubjects);
        
        %=== Extract info relative to each subject
        for j = 1: numSubjects
            r(j).Name = regexprep(s(j).subject, '\s+', ' '); % remove extra spaces in name and description
            r(j).Length = str2double(s(j).length);
            %r(j).Description = s(j).desc;
            
            subjectReport = queryReport(subjectStarts(j):subjectEnds(j));
            
            % hspsExpect = regexpi(subjectReport, 'Expect\(*\d*\)*\s+=\s{1,5}([e-\d\.]*)','tokens');
            hspsIdentities =  regexp(subjectReport, ...
                'Identities\s+=\s+(?<Match>(\d*))/(?<Possible>(\d*))\s+\((?<Percent>(\d*))','names');
            
            if use_positives
                hspsPositives =  regexp(subjectReport,...
                    'Positives\s+=\s+(?<Match>(\d*))/(?<Possible>(\d*))\s+\((?<Percent>\d*)%)','names');
            end
            if use_frame
                hspsFrames = regexp(subjectReport,'Frame\s*=\s*[+-]\d(\s*/\s*[+-]\d*)?','match');
            end
            if use_strand
                hspsStrands = regexp(subjectReport,...
                    'Strand\s*=\s*(Plus|Minus)\s*/\s*(Plus|Minus)','match');
            end
            
            [~, hspsScoreStarts] = regexpi(subjectReport, ...  removed: hspsScores
                'Score =\s{1,5}([e+\d\.]*)','tokens', 'start');
            numHSPs = length(hspsScoreStarts);
            
            if numHSPs > 1
                hspsScoreEnds = [hspsScoreStarts(2:end)-1 numel(subjectReport)];
            else
                hspsScoreEnds = numel(subjectReport);
            end
            
            %=== Preallocate space for each hsps
            dummy3.Score = [];
            dummy3.Expect = [];
            dummy3.Identities.Match = [];
            dummy3.Identities.Possible = [];
            dummy3.Identities.Percent = [];
            if use_frame
                dummy3.Frame = '';
            end
            if use_positives
                dummy3.Positives.Match = [];
                dummy3.Positives.Possible = [];
                dummy3.Positives.Percent = [];
            end
            if use_strand
                dummy3.Strand = '';
            end
            dummy3.Alignment = [];
            dummy3.QueryIndices = [];
            dummy3.SubjectIndices = [];
            
            rr = repmat(dummy3, 1, numHSPs);
            
            %=== Extract hsps info relative to subject j
            for k = 1:min(numHSPs,maxHSPs)
                % rr(k).Score = str2double(hspsScores{k}{:});
                % rr(k).Expect = str2double(regexprep(hspsExpect{k}{:}, '(^e[+-]\d)', '1$1')); % correct for abbreviated evalues
                rr(k).Identities.Match = str2double(hspsIdentities(k).Match);
                rr(k).Identities.Possible = str2double(hspsIdentities(k).Possible);
                rr(k).Identities.Percent = str2double(hspsIdentities(k).Percent);
                
                if use_positives
                    rr(k).Positives.Match = str2double(hspsPositives(k).Match);
                    rr(k).Positives.Possible = str2double(hspsPositives(k).Possible);
                    rr(k).Positives.Percent = str2double(hspsPositives(k).Percent);
                end
                if use_strand
                    rr(k).Strand = sscanf(char(hspsStrands{k}), 'Strand = %s %s %s');
                end
                if use_frame
                    rr(k).Frame = sscanf(char(hspsFrames{k}), 'Frame = %s %s %s');
                end
                
                %=== Extract alignment information (alignment, query and subject indices)
                hspsReport = subjectReport(hspsScoreStarts(k):hspsScoreEnds(k));
                alignments = regexp(hspsReport,['Query:*\s+(?<indq1>\d+)\s+(?<alnq>([a-zA-Z\-\*]+))\s+(?<indq2>\d+)\s', ...
                    '(?<alnc>([a-zA-Z\+|\s\*]+))\s+Sbjct:*\s+(?<inds1>\d+)\s+(?<alns>([a-zA-Z\-\*]+))\s+(?<inds2>\d+)'],'names');
                
                %=== Remove leading and trailing segments of \s that consensus might have
                for l = 1 : size(alignments,2)
                    actLen = numel(alignments(l).alnq);
                    consensus = char(regexprep(alignments(l).alnc, '\r$', ''));
                    lead = numel(consensus) - actLen + 1;
                    alignments(l).alnc = consensus(lead:lead + actLen - 1);
                end
                
                %=== Store relevant alignment info
                rr(k).Alignment =  [[alignments.alnq];[alignments.alnc];[alignments.alns]];
                rr(k).QueryIndices = [str2double(alignments(1).indq1), str2double(alignments(end).indq2)];
                rr(k).SubjectIndices = [str2double(alignments(1).inds1), str2double(alignments(end).inds2)];
                
            end % for k
            r(j).HSPs = rr;
            
        end % for j
        out(i).Hits = r;
        
    end % No hits found
    
    if numel(stats)==1
        out(i).Statistics = stats{1};
    else
        out(i).Statistics = stats{i};
    end
    times(i) = toc;
end % for i


%==========================================================================
function out = blastParser8(blasttext)
% Parse tabular BLAST reports (produced using local BLAST option m = 8 or 9)

if isempty(blasttext)
	warning(message('bioinfo:blastreadlocal:NoHits'));
	out = [];
	return
end

try
	info = textscan(blasttext,'%s %s %f %d %d %d %d %d %d %d %f %f');
	[~, reportEnds] = unique(info{1}, 'last');
	reportEnds = sort(reportEnds);    % undo sorting on queryNames, keep original order
	queryNames = info{1}(reportEnds);
	
	%=== Preallocate space for each query
	dummy1.Query = '';
	out = repmat(dummy1, 1, numel(reportEnds));
	
	
	%=== Parse report
	currLine = 1; % track lines as we read them
	reportStarts = [1; reportEnds + 1];
	for i = 1:numel(reportEnds) % for each query i
		out(i).Query = queryNames{i};
		[~, subjectEnds] = unique(info{2}(reportStarts(i):reportEnds(i)),'last');
		subjectEnds = sort(subjectEnds) + reportStarts(i) - 1;     % there is an offset for query i > 1
		subjectNames = info{2}(subjectEnds);
		numSubjects = numel(subjectNames);
		numHsps = diff([reportStarts(i)-1; subjectEnds]);
		
		%=== Preallocate space for each subject
		dummy2.Name = '';
		r = repmat(dummy2, 1,numSubjects);
		
		
		for j = 1 : numSubjects % for each subject j
			r(j).Name = subjectNames{j};
			
			%=== Preallocate space for each hsp
			dummy3.AlignmentLength = [];
			dummy3.Score = [];
			dummy3.Expect = [];
			dummy3.Identities.Match = [];
			dummy3.Identities.Possible = [];
			dummy3.Identities.Percent = [];
			dummy3.Mismatches.Match = [];
			dummy3.Mismatches.Possible = [];
			dummy3.Mismatches.Percent = [];
			dummy3.Gaps.Match = [];
			dummy3.Gaps.Possible = [];
			dummy3.Gaps.Percent = [];
			dummy3.QueryIndices = [];
			dummy3.SubjectIndices = [];
			rr = repmat(dummy3, 1, numHsps(j));
			
			
			for k = 1:numHsps(j)
				rr(k).AlignmentLength = info{4}(currLine);
				rr(k).Score = info{12}(currLine);
				rr(k).Expect = info{11}(currLine);
				rr(k).Identities.Match = info{3}(currLine) * rr(k).AlignmentLength / 100;
				rr(k).Identities.Possible = rr(k).AlignmentLength;
				rr(k).Identities.Percent = info{3}(currLine);
				rr(k).Mismatches.Match = info{5}(currLine);
				rr(k).Mismatches.Possible = rr(k).AlignmentLength;
				rr(k).Mismatches.Percent = rr(k).Mismatches.Match * 100 /rr(k).AlignmentLength;
				rr(k).Gaps.Match = info{6}(currLine);
				rr(k).Gaps.Possible = rr(k).AlignmentLength;
				rr(k).Gaps.Percent = rr(k).Gaps.Match * 100 /rr(k).AlignmentLength;
				rr(k).QueryIndices = [info{7}(currLine) info{8}(currLine)];
				rr(k).SubjectIndices = [info{9}(currLine) info{10}(currLine)];
				currLine = currLine + 1;
			end % for k
			r(j).HSPs = rr;
		end % for j
		out(i).Hits = r;
	end % for i
catch lerr
	error(message('bioinfo:blastreadlocal:CorruptedFormat'));
end
