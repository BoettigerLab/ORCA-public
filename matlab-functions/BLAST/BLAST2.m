function [blastData] = BLAST2(sequence,dataBase,varargin)
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
defaults(end+1,:) = {'maxhits', 'positive', 30};
defaults(end+1,:) = {'maxResults', 'positive', 30};
defaults(end+1,:) = {'showplots', 'boolean', false};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'numThreads', 'positive', 1};

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
        parameters.blastPath = 'C:\"Program Files"\NCBI\LegacyBLAST\bin\';
    else
        parameters.blastPath = 'C:\"Program Files"\NCBI\blast-2.2.27+\bin\';
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
if ~isstruct(sequence)
    fastafile = sequence;
else
    fastafile = [scratchPath,'tempfasta.fasta'];
    if exist(fastafile)
        delete(fastafile)
    end
    WriteFasta(fastafile,sequence);
end

%--------------------------------------------------------------------------
% Overwrite output file if it already exists
%--------------------------------------------------------------------------
if exist(parameters.outputFile,'file') == 2
    disp(['deleting ',parameters.outputFile]);
    del = system(['del ', parameters.outputFile]);
    dd = 0;
    while exist(parameters.outputFile,'file') && ~del
        disp('file in use, making new temp file...');
        dd = dd + 1;
        parameters.outputFile = regexprep(parameters.outputFile,'temp',['temp',num2str(dd)]);
        del = system(['del ', parameters.outputFile]);
    end
end

%--------------------------------------------------------------------------
% Run blast
%--------------------------------------------------------------------------
if ~parameters.legacy
    command = [parameters.blastPath 'blastn.exe' ...
        ' -query ' fastafile ...
        ' -task ' ' "blastn-short" ' ...
        ' -db ' dataBase ...
        ' -out ' parameters.outputFile ...
        ' -num_alignments ' num2str(parameters.maxResults) ...
        ' -num_threads ' num2str(parameters.numThreads) ...
        ' -outfmt ' num2str(0)];
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' command]);
    display('-----------------------------------------------------------------');
    system(command);
else    
    command = [parameters.blastPath 'blastall.exe' ...
        ' -i ' fastafile ...
        ' -p ' 'blastn' ...
        ' -d ' dataBase ...
        ' -o ' parameters.outputFile ...
        ' -K ' num2str(parameters.maxResults)];
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' command]);
    display('-----------------------------------------------------------------');
    system(command);
end

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
blastData = blastreadlocal(parameters.outputFile, 0);

    
    


