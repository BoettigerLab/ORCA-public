function [probeFasta,targetRegions,probeSets] = BuildChrTracingLibRNA(analysisSavePath,locus,stepSize,locusGeneName,varargin)
% -------------------------------------------------------------------------
% Function help file:
% -------------------------------------------------------------------------
% BuildChrTracingLib(requiredPar1,requiredPar2) 
%
% Description of what my function does when passed requiredPar1 and
% requiredPar2.  
%
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% requiredPar1, datatype, description of what this par is 
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% output1, datatype, description of what this par is 
% 
% -------------------------------------------------------------------------
% Optional Inputs
% -------------------------------------------------------------------------
% format: 'name', 'datatype', defaultValue, description 
%
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
% 
% Alistair Boettiger
% Nov 17, 2016
% 
% END of help file.


% -------------------------------------------------------------------------
% Defaults for optional variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; % example optional par and its default 
defaults(end+1,:) = {'fwdPrimers', 'fasta', fastaread('C:\Data\Oligos\FwdPrimers.fasta')};
defaults(end+1,:) = {'revPrimers', 'fasta', fastaread('C:\Data\Oligos\RevPrimers.fasta')};
defaults(end+1,:) = {'secondaries', 'fasta', fastaread('C:\Data\Oligos\ReadoutSeqs.fasta')};
defaults(end+1,:) = {'commonRT', 'string', 'catcaacgccacgatcagct'};  % 20 bp complimentary to P4-405 end. 
defaults(end+1,:) = {'genomeFasta', 'string', 'C:\Data\Fly\dm3\dm3.fasta'}; 
defaults(end+1,:) = {'flyRepeatFasta', 'string', 'C:\Data\Fly\Dmel_repeats.fasta'}; 
defaults(end+1,:) = {'genomeOffTarget','boolean',true};
defaults(end+1,:) = {'startPrimer','integer',1};
defaults(end+1,:) = {'offsetPrimerPairs','integer',1};
defaults(end+1,:) = {'minProbesPerRegion','integer',15};
defaults(end+1,:) = {'probeLength','integer',40};
defaults(end+1,:) = {'TmRange','array',[65, 90]};
defaults(end+1,:) = {'reverseTile','boolean',false};
defaults(end+1,:) = {'parallel','integer',4};
defaults(end+1,:) = {'repeatOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'secOTrange','array', [-1 0]};  % [-1 0])
defaults(end+1,:) = {'repeatOTmatch','integer', 14};  % [-1 0])
defaults(end+1,:) = {'secOTmatch','integer', 12};  % [-1 0])
% defaults(end+1,:) = {'optionalVar', '<VarType>', <VarValue>};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'a requiredPar1 (datatype) and requiredPar2 (datatype) are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
%% Actual Function
%--------------------------------------------------------------------------

if parameters.verbose
    disp('creating probes for:')
    disp(locus)
end

% a temporary fasta file which will be overwritten by successive function calls 
regionFasta = [analysisSavePath,'RNA_Regions_temp.fasta'];

% need to get the fasta-headers in the correct format
[rnaNames,rnaSeqs] =  GetTranscribedSeqs(locus) ;
for r=1:length(rnaNames)
    rnaNames{r} = [locusGeneName,'_',locus,' gene=',rnaNames{r}];
end
try
    WriteFasta(regionFasta,rnaNames,rnaSeqs,'Append',false,'Warnings',true);
catch er
    disp(er.message);
    disp(rnaNames);
    disp(rnaSeqs);
    return
end

% Build off-target fasta from the rest of the genome.
% (consider just using the chromosome containing the region of interest)
if parameters.genomeOffTarget
    regions = {locus}; 
    otFasta = BuildOffTargetFasta(regions,parameters.genomeFasta);
    WriteFasta([analysisSavePath,'offTarget.fasta'],otFasta);
end

[probeFasta,targetRegions,probeSets] = FastaToSeqProbes(regionFasta,analysisSavePath,'parameters',parameters,'locusGeneName',locusGeneName); 



