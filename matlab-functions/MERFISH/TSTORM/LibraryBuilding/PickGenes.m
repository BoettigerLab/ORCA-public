function pickedGenes = PickGenes(ProbeData,numGenes,numOligos,varargin)
%--------------------------------------------------------------------------
% Pick genes randomly based on user requirements
%
% Required Parameters
% ProbeData - structure
% 
%% DefaultParameters
% RequiredGenes = {};  
% BinnedFPKM = 0;  number of bins between min and max FPKM
% minFPKM = 0;


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'requiredGenes', 'cell',{} };
defaults(end+1,:) = {'excludedGenes', 'cell',{} };
defaults(end+1,:) = {'primerNum', 'cell',{} };
defaults(end+1,:) = {'desiredPerBin', 'positive',[] };
defaults(end+1,:) = {'binsFPKM','positive',1};
defaults(end+1,:) = {'minFPKM', 'nonnegative',0 };
defaults(end+1,:) = {'maxFPKM', 'positive',inf };
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 3
    error('matlabSTORM:invalidArguments', 'required: ProbeData,numGenes,numOligos');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

requiredGenes = parameters.requiredGenes;
excludedGenes = parameters.excludedGenes;
desiredPerBin = parameters.desiredPerBin;
binsFPKM =      parameters.binsFPKM;
minFPKM =       parameters.minFPKM;
maxFPKM =       parameters.maxFPKM; 

%-------------------------------------------------------------------------
%% Main Function
%-------------------------------------------------------------------------
if ~isempty(desiredPerBin) && length(desiredPerBin) < binsFPKM
    error('length of desiredPerBin must equal value of binnedFPKM');
end

% A little repair to the ProbeData structure
noProbes = cellfun(@isempty, ProbeData.Nprobes) ;
ProbeData.Nprobes(noProbes) = {0}; 
ProbeData.CommonName(noProbes) = '';

% cut on min and max FPKM
goodFPKM = ([ProbeData.FPKM(:)] < maxFPKM & [ProbeData.FPKM(:)] > minFPKM)';

% cut on number of probes
sufficientProbes = [ProbeData.Nprobes{:}] > numOligos;

% Find indices of Required Genes
reqIdx = StringFind(ProbeData.CommonName,requiredGenes,'exactly',true)';
if iscell(reqIdx)
    IsoNames = cellfun(@(x) x(1:end-2),ProbeData.IsoformName,'UniformOutput',false);
    reqIdx = StringFind(IsoNames,requiredGenes,'exactly',true)';
end

% Find indices of Excluded Genes
exIdx = zeros(length(excludedGenes),1); 
for g=1:size(excludedGenes,1)
    exIdx(g) = find(strcmp(cellstr(char(ProbeData.CommonName{:})),excludedGenes{g}));
end

% Indices of all genes matching search characteristics
viableGenes = goodFPKM & sufficientProbes;
viableGenes(ismember(viableGenes,exIdx)) = []; % remove excluded genes

deficit = 0;
if binsFPKM > 0
    rangeFPKM = logspace(log10(minFPKM),log10(maxFPKM),binsFPKM+1);
    numPerBin = floor((numGenes)/binsFPKM);
    idxSelect = cell(binsFPKM,1); 
    for n=1:binsFPKM % n = 3
        minLevel =  rangeFPKM(binsFPKM-n+1);
        maxLevel = rangeFPKM(binsFPKM-n+2);
        inRange = find([ProbeData.FPKM(:)] > minLevel & [ProbeData.FPKM(:)] <= maxLevel & [ProbeData.Nprobes{:}]' > numOligos );
        inRange(ismember(inRange,exIdx)) = []; % Remove exluded genes
        inRange = inRange(randperm(length(inRange))); % randomly permute indices
        numFound = length(inRange);
        numChosen = min(numPerBin+deficit,numFound);
        
        if ~isempty(desiredPerBin)
            numPerBin = desiredPerBin(n);
            numChosen = min(numChosen,numPerBin);
        end

        % If we don't find enough in this bin, grab extra genes from the next
        % bin to make up the deficit.  
        if numFound < numPerBin + deficit
            warning(['Found only ',num2str(numFound),' genes with FPKMS ',num2str(minLevel),'-',num2str(maxLevel),'.  Desired ',num2str(numPerBin+deficit),'. Will use genes from next bin to fill the deficit']); 
            deficit = deficit + (numPerBin-numFound);
        else
            deficit = 0;
        end
        idxSelect{n} = inRange(1:numChosen);
        if isempty(idxSelect{n})
            idxSelect{n} = [];
        end       
    end
end

 pickedIdx = cat(1,idxSelect{:});

 
% if required genes aren't by some chance already included, add them to the
% top
reqIdxUsed = ismember(reqIdx,pickedIdx); %
if ~isempty( reqIdx(~reqIdxUsed) )
    pickedIdx = [reqIdx(~reqIdxUsed); pickedIdx];
end

% if we have too many genes truncate
if length(pickedIdx) > numGenes
    pickedIdx = pickedIdx(1:numGenes); 
end

% Fill out any extra space with randomly selected viable Genes
if length(pickedIdx) < numGenes
    disp(['Not enough genes available from desired bins. Adding ',num2str(numGenes-length(pickedIdx)),' more genes at random']); 
    newIdx = setdiff(find(viableGenes),pickedIdx);
    newIdx = newIdx(randperm(length(newIdx)));
    newIdx = newIdx(1:numGenes - length(pickedIdx))';
    pickedIdx = cat(1,newIdx,pickedIdx);
end
    
pickedGenes = IndexStructure(ProbeData,pickedIdx,'celldata',false);
% newPickedGenes    = IndexStructure(ProbeData,pickedIdx,'celldata',false);
    
