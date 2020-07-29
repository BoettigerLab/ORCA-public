function [polymer,distMap] = TableToPolymer(dTable,varargin)
% 
% Update History
% adapted from TableToDistMap to create a sleaker, non-backwards compatible
% form of the function. 

defaults = cell(0,3);
defaults(end+1,:) = {'dims',{2,3},3};
defaults(end+1,:) = {'nSteps','freeType',[]};
defaults(end+1,:) = {'chromCorrect','boolean',false};
defaults(end+1,:) = {'sortMethod',{'byHybe','byRead'},'byHybe'};
defaults(end+1,:) = {'dataBasedDriftCorrect','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% save some speed
if nargout>1 
    computeDistMap = true;
else
    computeDistMap = false;
end
% avoid matrix concatination trouble by allowing the total number of steps
% in the matrix to be passed as a variable


spotIDs = unique(dTable.s); % as a row vector
nPolys = length(spotIDs);

%% Main function
%
switch pars.sortMethod
    case 'byHybe'
        if isempty(pars.nSteps)
            numHybes=max(dTable.hybe);
        else
           numHybes = pars.nSteps;
        end
        % initialize the matrices
        polymer = nan(numHybes,pars.dims,nPolys);
        distMap = nan(numHybes,numHybes,nPolys);       
        % populate matrices and compute distances
        for s=1:nPolys
            sTable = dTable(dTable.s==spotIDs(s),:);
            if ~isempty(sTable)
                    polymer(sTable.hybe,:,s) = [sTable.x,sTable.y,sTable.z];  
                if computeDistMap
                    distMap(:,:,s) = squareform(pdist(polymer(:,1:pars.dims,s)));
                end
            end
        end
    case 'byRead'
        % organize by readout number
        %   computing the controls is a bit annoying this way needs work
        %   later. 
        isDat = strcmp(dTable.dataType,'H');
        datTable = dTable(isDat,:);
        if isempty(pars.nSteps)
            nReads = max.datTable.readout;
        else
           nReads = pars.nSteps;
        end
        % initialize the matrices
        polymer = nan(nReads,pars.dims,nPolys);
        distMap = nan(nReads,nReads,nPolys);
        % populate matrices and compute distances
        for s=1:nPolys    
            sTable = datTable(datTable.s==spotIDs(s),:);
            if ~isempty(sTable)
                polymer(sTable.readout,:,s) = [sTable.x,sTable.y,sTable.z]; 
                if computeDistMap
                    distMap(:,:,s) = squareform(pdist(polymer(:,1:pars.dims,s)));
                end
            end
        end
    otherwise
        error(['sort method ',pars.sortMethod,' not recognized. Valid values: "byHybe" or "byRead"']);
end

if pars.chromCorrect
   error('not implemented here yet');  
end

