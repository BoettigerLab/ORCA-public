% Demo PrimerDesigner
% 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% License: 
% Demo by Alistair Boettiger
% Built on the PrimerDesigner class written by Jeffrey Moffitt, May 4, 2015
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% restoredefaultpath

% First we need to add 'OligoBuilding' folder to our file paths. 
% If we run this script, Matlab will change the working directory to the 
% DemoScripts folder.  We then ask matlab to return a list of 
matlab_public_functions = 'C:\Users\Alistair\Desktop\code\matlab-public-functions'; % UPDATE THIS
saveDataFolder = 'C:\Data';  % UPDATE THIS
% addpath(genpath(matlab_public_functions)); 

% Path to the fasta files which include the sequences we want to avoid
offTargetSeqs =[matlab_public_functions,'\DemoData\OffTargetSeqs.fasta'];

% Add an off target table
OT15 = OTTable(fastaread(offTargetSeqs), 15, 'verbose', true);

% Add an additional off-target table with a different homology cut-off
OT17 = OTTable(fastaread(offTargetSeqs), 17, 'verbose', true);

% Design Primers
prDesigner = PrimerDesigner('numPrimersToGenerate', 10, ...
        'OTTables', [OT15,OT17], ...  % pass the tables
        'OTTableNames', {'AbundGenes','LessAbundGenes'},...
        'ntComposition',[0.25 0.25 0.25 0.25]); % assign some names to the off-target tables for book-keeping
   
int2nt(prDesigner.seqs+1) % convert integer seqs to ATGC (note the +1) 
    
% We could pass other parameters such as
%   {'parameterFlag', 'dataType', default value}
%   {'verbose', 'boolean', true}; % Display progress of construction
%   {'ntComposition', 'positive', [0.25 0.25 0.25 0.25]};
%   {'OTTables', 'freeType', []};
%   {'OTTableNames', 'cell', {}};
%   {'parallel', 'parallel', []};
%   {'seqs', 'array', []};
%   {'primerLength', 'positive', 20};
%   {'numPrimersToGenerate', 'positive', 1e6};
%   {'homologyMax', 'positive', 8};
%   {'monovalentSalt', 'positive', 0.3};
%   {'primerConc', 'positive', 0.5e-6};
%   {'seqsToRemove', 'cell', {'AAAA', 'TTTT', 'GGGG', 'CCCC'}};
    
prDesigner.WriteFasta([saveDataFolder,'\testPrimers.fasta'],'namePrefix','primer');