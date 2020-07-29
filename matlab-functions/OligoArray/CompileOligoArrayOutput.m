function CompileOligoArrayOutput(oligo_folder,varargin)
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% oligo_folder - string. folder containing the Output files of OligoArray
% should conatin the endings _oligos.txt.
% 
% -------------------------------------------------------------------------
% Optional Inputs
% -------------------------------------------------------------------------
% 'savePath' folder to save ProbeData.mat in.
% 'blastLib' list of sequences BLASTED against

% exampleData = 'C:\Users\Alistair\Documents\Research\Projects\MERFISH-Release\MERFISH-data\MERFISH_Examples\';
% dataFolder = [exampleData,'probe_construction_data\'];
% oligo_folder = [exampleData,'probe_construction_output\'];
% blastLib =  [dataFolder,'TargetGeneSeqs.fasta'];

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'savePath', 'string', ''};
defaults(end+1,:) = {'blastLib', 'string', ''};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'required: oligo_folder');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


if isempty(parameters.savePath)
    savePath = oligo_folder;
else
    savePath = parameters.savePath;
end

%% 

% Parsing variable inputs
if ~isempty(parameters.blastLib)
    blastSeqs = fastaread(parameters.blastLib);
    allNames = {blastSeqs.Header}';
else
    allNames = {''};
end

% Scan for files and load them
disp('Searching for _oligo.txt files...'); 
oligofiles = dir([oligo_folder,'*_oligos.txt']);
numGenes = length(oligofiles);
disp(['Found oligo sets for ', num2str(numGenes), ' genes']);

% make sure files completed running and saved a log file.
haslog =  dir([oligo_folder,'*_log.txt']);
missinglog = numGenes - length(haslog);
disp(['OligoArray failed for ',num2str(missinglog),' genes']);
    
clear ProbeData;
numGenes =length(oligofiles) ; % 100
disp('screening oligos... ');   

% initialize array
ProbeData.Nprobes = cell(numGenes,1);
ProbeData.GeneName = cell(numGenes,1);
ProbeData.FivePrimeEnd  = cell(numGenes,1);
ProbeData.FreeEnergy  = cell(numGenes,1);
ProbeData.Enthalpy  = cell(numGenes,1);
ProbeData.Entropy  = cell(numGenes,1);
ProbeData.MeltingTemp  = cell(numGenes,1); 
ProbeData.Sequence = cell(numGenes,1);
    
for n = 1:numGenes  % n = 1  % n= 4851
    disp(['processing gene ',num2str(n),' of ',num2str(numGenes)]);
    oligofileName = oligofiles(n).name;
    oligofile = [oligo_folder,oligofileName];
    
    fid = fopen(oligofile);
    fmt = '%s %f %f %f %f %f %f %s %s %s %*[^\n]';
    data = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','TreatAsEmpty','NA');
    fclose(fid);

    if ~isempty(data{1})
        % Let's get target gene name from the header
        genename = regexprep(oligofileName,'_oligos.txt','');
        
        %  Search for off-target hits recorded in the file
        Nprobes = length(data{3}(:,2));
        offTargets = false(Nprobes,1); 
        currTargetID = StringFind(allNames,genename,'boolean',true);
        offTargetNames = allNames(~currTargetID);
        for j=1:Nprobes  % j = Nprobes  j =2
             hits = StringFind(data{3}(j,1),offTargetNames,'cellOutput',true);
             offTargetHits = cat(1,hits{:});
             if ~isempty(offTargetHits)
                 offTargets(j) = true;
             end
        end
        
        % Record probe properties into a structure array
        ProbeData.Nprobes{n} = length(data{3}(~offTargets,2));
        ProbeData.GeneName{n} = genename;
        ProbeData.FivePrimeEnd{n} = data{2}(~offTargets,1);
        ProbeData.FreeEnergy{n} = data{2}(~offTargets,3);
        ProbeData.Enthalpy{n} = data{2}(~offTargets,4);
        ProbeData.Entropy{n} = data{2}(~offTargets,5);
        ProbeData.MeltingTemp{n} =data{2}(~offTargets,6);  
        ProbeData.Sequence{n} = data{3}(~offTargets,2);

    else
        % add appropriate empty delimiters if data is missing
        ProbeData.Nprobes{n} =0;
        ProbeData.GeneName{n} = genename;
        ProbeData.FivePrimeEnd{n} = NaN;
        ProbeData.FreeEnergy{n} = NaN;
        ProbeData.Enthalpy{n} = NaN;
        ProbeData.Entropy{n} = NaN;
        ProbeData.MeltingTemp{n} =NaN;
        ProbeData.Sequence{n} = '';
    end
end      

save([savePath,'ProbeData.mat'],'ProbeData'); 
disp(['wrote ',savePath,'ProbeData.mat']);

