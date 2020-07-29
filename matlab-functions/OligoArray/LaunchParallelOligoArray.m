function BatchLaunchOligoArray(oligoArrayCommandtemp,geneFasta,varargin)
% takes an oligoArray command

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'savePath', 'string', ''};
defaults(end+1,:) = {'headerSeps', 'string',''};
defaults(end+1,:) = {'maxFragment', 'positive',1E3};
defaults(end+1,:) = {'batchsize', 'positive',4 };
defaults(end+1,:) = {'verbose', 'boolean',true };
defaults(end+1,:) = {'preprun','boolean',false};

if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires an oligoArrayCommand string and geneGasta');
end

parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% parsing input a bit more
if isempty(parameters.savePath) 
    s = strfind(oligoArrayCommandtemp,'-o ');
    e = strfind(oligoArrayCommandtemp,'-r ');
    parameters.savePath = extractpath(oligoArrayCommandtemp(s+3:e-3)); 
end


%% Launch jobs

% First split fasta files into chunks and force gene names into OligoArray
% compatable names. 
Gs = length(geneFasta);
isoNames = cell(Gs,1);
for g=1:Gs
    locusName = geneFasta(g).Header;
    sequence = geneFasta(g).Sequence;
    
    % Trim gene_name and remove illegal characters
    locusName = regexprep(locusName,...
        {'\.','\:','\s',',','\\','\/'},{'P','C','_','',']','['}); 
    if ~isempty(parameters.headerSeps)
        sps = strfind(geneFasta(g).Header,parameters.headerSeps);
        locusName = locusName(sps(1)+length(parameters.headerSeps):sps(2)-1);
    end
        
    % Decompose region into smaller fragments
    seq_length = length(sequence); 
    locusSaveName =  locusName;
    gene_fragment_seqs = cell(100,1);
    gene_fragment_names = cell(100,1); 
   subfragment_start = 1;
   subfragment_end = min(seq_length,parameters.maxFragment);
   f = 0; % fragment counter
   while  subfragment_start < seq_length 
        f = f+1;
        gene_fragment_seqs{f} = ...
        sequence(subfragment_start:subfragment_end);
        subfragment_start = subfragment_start + parameters.maxFragment;
        subfragment_end = min(subfragment_end + parameters.maxFragment,seq_length);
        gene_fragment_names{f} = [locusName,'_pt',num2str(f)];
   end

    hasdata = logical(1-cellfun(@isempty,gene_fragment_seqs));
    gene_fragment_names = gene_fragment_names(hasdata);
    gene_fragment_seqs = gene_fragment_seqs(hasdata);
    Gene_fasta.Sequence = gene_fragment_seqs;
    Gene_fasta.Header = gene_fragment_names;

    fileout = [parameters.savePath,locusSaveName,'.fasta'];
    WriteFasta(fileout,Gene_fasta,[],'Append',false,'Warnings',false); 
    isoNames{g} = locusSaveName;
end  

verbose = parameters.verbose; 

if parameters.preprun
    oligoArrayCommands = cell(Gs,1); 
    for g=1:Gs
        locusSaveName = isoNames{g};
        oligoArrayCommands{g} = regexprep(oligoArrayCommandtemp,'genename',locusSaveName);
    end
    save([parameters.savePath,'oligoArrayCommands.mat'],'oligoArrayCommands');
    if verbose
        disp(['wrote ',parameters.savePath,'oligoArrayCommands.mat']);
    end
    
else
    parpool(parameters.batchsize)
    parfor g=1:Gs
        % Replace 'fastainput' with the fasta name for this gene.
        % Then send OligoArrayCommand to Queue. 
        locusSaveName = isoNames{g}
        oligoArrayCommand = regexprep(oligoArrayCommandtemp,'genename',locusSaveName);
        system(oligoArrayCommand);

        if verbose
         display('-----------------------------------------------------------------');
         display(['Running file ',num2str(g),' of ',num2str(Gs),':'])
         display(['     ' oligoArrayCommand]);
         display('-----------------------------------------------------------------');
        end 
    end
end
