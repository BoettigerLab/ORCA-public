function [blastResults,blastData,allHits] = BLAST(sequence,dataBase,varargin)
% ------------------------------------------------------------------------
% [blastResults,blastData,allhits] = BLAST(sequence,dataBase,varargin)
% This function blasts the fastaFile specified by sequence against the
% dataBase.
%--------------------------------------------------------------------------
% Necessary Inputs
% fastaFile: String to a fasta file or a fasta structure
% dataBase: Path to a valid database
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global variables
%-------------------------------------------------------------------------- 
global scratchPath

%--------------------------------------------------------------------------
% Default Parameters
%-------------------------------------------------------------------------- 
defaults = cell(0,3);
defaults(end+1,:) = {'blastPath', 'string', ''}; % Base tag for all images
defaults(end+1,:) = {'legacy', 'boolean', false};
defaults(end+1,:) = {'outputFile', 'string', [scratchPath '\blastScratch.txt']};
defaults(end+1,:) = {'exludeFirstHit', 'boolean', false};
defaults(end+1,:) = {'no3pMatch', 'boolean', false};
defaults(end+1,:) = {'primers', 'struct', []};
defaults(end+1,:) = {'maxhits', 'positive', 100};
defaults(end+1,:) = {'maxResults', 'positive', 100};
defaults(end+1,:) = {'wordSize','positive',15}; 
defaults(end+1,:) = {'showplots', 'boolean', false};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'numThreads','integer',1};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabFunctions:invalidArguments', 'A valid fasta file (or structure) is required and a database');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%-------------------------------------------------------------------------
% Handle default blast paths: Backwards compatibility for Alistair
%-------------------------------------------------------------------------
if isempty(parameters.blastPath)
    if parameters.legacy
        parameters.blastPath = 'C:\Software\NCBI\LegacyBLAST\bin\';
    else
        parameters.blastPath = 'C:\Software\NCBI\blast-2.15.0+\bin\';
    end
end

%--------------------------------------------------------------------------
% Test Installation
%--------------------------------------------------------------------------
% Test BLAST installations
if ~parameters.legacy
    if ~exist([regexprep(parameters.blastPath,'"','') 'makeblastdb.exe'],'file')
        error(['Could not find NCBI blast+.  Please update blastPlusPath in ',mfilename]);
    end
else
    if ~exist([regexprep(parameters.blastPath,'"','') 'formatdb.exe'],'file')
        error(['Could not find NCBI blast.  Please update legacyBlastPath in ',mfilename]);
    end
end

%--------------------------------------------------------------------------
% Write Fasta file if needed
%--------------------------------------------------------------------------
% if a filepath is passed
if ~isstruct(sequence) && ischar(sequence) && exist(sequence,'file') == 2
    fastafile = sequence; 
    writeFasta = false;
elseif ~isstruct(sequence) && ischar(sequence)  % if a sequence is passed  
    fa.Header = 'InputSequence';
    fa.Sequence = sequence;
    writeFasta = true;
elseif isstruct(sequence) % if a matlab fasta structure is past
    fa = sequence;
    writeFasta = true;
end
if writeFasta % write a fasta file if necessary
    fastafile = [scratchPath,'tempfasta.fasta'];
    if exist(fastafile,'file')~=0
        [~,cmdOut] = system(['del ',fastafile]);
        if ~isempty(cmdOut) % contains(cmdOut(2,:),'cannot access')
            if parameters.verbose
                warning('failed to delete existing tempfasta, creating new tempfasta');
            end
            fastafile = IncrementSaveName(fastafile);
        else
            if parameters.verbose
                disp('rewriting temp fasta file...')
            end
        end
    end
    fastawrite(fastafile,fa);
    %WriteFasta(fastafile,fa,[]);
end

%--------------------------------------------------------------------------
% Overwrite output file if it already exists
%--------------------------------------------------------------------------
if exist(parameters.outputFile,'file') == 2
    if parameters.verbose
        disp(['deleting ',parameters.outputFile]);
    end
    [~,cmdOut] = system(['del ', parameters.outputFile]);
    del = isempty(cmdOut);
    dd = 0;
    while exist(parameters.outputFile,'file') && ~del
        if parameters.verbose
            disp('file in use, making new temp file...');
        end
        dd = dd + 1;
        parameters.outputFile = regexprep(parameters.outputFile,'temp',['temp',num2str(dd)]);
        [~,cmdOut] = system(['del ', parameters.outputFile]);
        del = isempty(cmdOut);
    end
end

%--------------------------------------------------------------------------
% Run blast
%--------------------------------------------------------------------------
if ~parameters.legacy
    dataBase = regexprep(dataBase,'.fasta','');
    command = [parameters.blastPath 'blastn.exe' ...
        ' -query ' fastafile ...
        ' -task ' ' "blastn-short" ' ...
        ' -db ' dataBase ...
        ' -out ' parameters.outputFile ...
        ' -num_alignments ' num2str(parameters.maxResults) ...
        ' -num_threads ' num2str(parameters.numThreads) ...
        ' -word_size ',num2str(parameters.wordSize) ...
        ' -outfmt ' num2str(0) ...
        ' &'];
    if parameters.verbose
        disp('-----------------------------------------------------------------');
        disp('Issuing:')
        disp(['     ' command]);
        disp('-----------------------------------------------------------------');
    end
    % system(command);
    
    % Run command and wait for completion
    prc = SystemRun(command,'Hidden',false); 
    runFinished = prc.HasExited;
    while ~runFinished
        pause(1);
        runFinished = prc.HasExited;
    end
    
else    
    command = [parameters.blastPath 'blastall.exe' ...
        ' -i ' fastafile ...
        ' -p ' 'blastn' ...
        ' -d ' dataBase ...
        ' -o ' parameters.outputFile ...
        ' -a ' num2str(parameters.numThreads) ...
        ' -K ' num2str(parameters.maxResults)];
    if parameters.verbose
        disp('-----------------------------------------------------------------');
        disp('Issuing:')
        disp(['     ' command]);
        disp('-----------------------------------------------------------------');
    end
    system(command);
end

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
if parameters.verbose
    disp('loading data...');
end

blastData = BLASTreadLocal(parameters.outputFile, 0); % 0  % My version maxes out on mapping super repetitive regions  
% blastData = blastreadlocal(parameters.outputFile, 0); % 0



% % DEPRECIATED.  See Hits per Query below
% if parameters.maxhits > 0
%     [allhits] = MaxBLASThits(blastData,'verbose',false,...
%         'no3pMatch',parameters.no3pMatch,...
%         'primers',parameters.primers,...
%         'maxhits',parameters.maxhits,...
%         'exludeFirstHit',parameters.exludeFirstHit);
%     if parameters.showplots
%         hist(allhits);
%         xlabel('length of homology region'); 
%     end
% else
%     allhits = [];
% end




blastResults = ParseBlastData(blastData); 

if parameters.verbose && ~isempty(blastData(1).Hits)
    fprintf('\n'); 
    disp('top BLAST hit per query:')
    disp(  [blastResults.topHitName,blastResults.topHitSeq] );
end


try
    if ~isempty(blastData(1).Hits)
        topHit(1).Name = blastData.Hits(1).Name;
        topHit(1).Sequence = blastData.Hits(1).HSPs(1).Alignment;
        if parameters.verbose
            fprintf('\n'); 
            disp('top BLAST hit:')
            disp(topHit(1).Name);
            disp(topHit(1).Sequence);
        end
        if length(blastData.Hits) > 1 && parameters.maxhits > 1
            topHit(2).Name = blastData.Hits(2).Name;
            topHit(2).Sequence = blastData.Hits(2).HSPs(1).Alignment;
            if parameters.verbose
                fprintf('\n'); 
                disp('2nd BLAST hit:')
                disp(topHit(2).Name);
                disp(topHit(2).Sequence);
            end      
        end
    end   
catch er
    if parameters.verbose
        disp(er.getReport)
    end
end


% Hits per querry
numProbes = length(blastResults.numHits);
allHits = zeros(numProbes,1);
hitLocs = cell(numProbes,1);
for p=1:numProbes
    totHits = 0;
    numChrsHit = length(blastData(p).Hits);
    hitLocs{p} = cell(numChrsHit,1);
    for c=1:numChrsHit
        totHits = totHits + length(blastData(p).Hits(c).HSPs);
        hitLocs{p}{c} = cat(1,blastData(p).Hits(c).HSPs.SubjectIndices);
    end
    allHits(p) = totHits;
end   
    


