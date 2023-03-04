function [polymer,distMap,spotData] = TableToPolymer(dTable,varargin)
% Inputs
% dTable - either the path to an AllFits table or the loaded table.
% 
% Outputs 
% polymer, distMap, spotData
%  
% Update History
% adapted from TableToDistMap to create a sleaker, non-backwards compatible
% form of the function. 
%
% Note, it is faster to parse data by fov and then recombine. 

defaults = cell(0,3);
defaults(end+1,:) = {'dims',{'xyz','xy','xz'},'xyz'};
defaults(end+1,:) = {'bins','nonnegative',0}; % 0 for auto detect
defaults(end+1,:) = {'chromCorrect','boolean',true};
defaults(end+1,:) = {'computeDistMap','boolean',true};
defaults(end+1,:) = {'sortMethod',{'byHybe','byRead'},'byRead'};
defaults(end+1,:) = {'selectChn','nonnegative',0}; % 0 for all chn
defaults(end+1,:) = {'dataBasedDriftCorrect','boolean',false};
defaults(end+1,:) = {'removeBlank','boolean',false};
defaults(end+1,:) = {'shiftSign','integer',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% save some speed
computeDistMap = pars.computeDistMap;
if nargout<2 
    computeDistMap = false;
end

if ischar(dTable) || isstring(dTable)
    dTable = readtable(dTable);
end

if ~strcmp(dTable.Properties.VariableNames,'xcShift')
    pars.chromCorrect = false;
end

% load table if value is a string
if pars.selectChn ~= 0
    dTable = dTable(dTable.chn == pars.selectChn,:);
end

% avoid matrix concatination trouble by allowing the total number of steps
% in the matrix to be passed as a variable
distMap = [];
spotData = [];
spotIDs = unique(dTable.fs); % as a row vector
nPolys = length(spotIDs);

%% Main function
%

switch pars.dims
    case 'xyz'
        dims = 1:3;
    case 'xy'
        dims = 1:2;
    case 'xz'
        dims = [1,3];
end


% data based drift correction
if pars.dataBasedDriftCorrect
    % remove fiduical based drift correction
    if quantile(abs(dTable.xshift),.9) < 100
        % old version in pixels
       dTable.x = dTable.x + dTable.xshift*154;  % this was an old table format with shifts recorded in pixels but x,y,z recorded in nm. Now everything is nm 
       dTable.y = dTable.y + dTable.yshift*154;
       dTable.z = dTable.z + dTable.zshift; 
    else
       dTable.x = dTable.x + dTable.xshift;
       dTable.y = dTable.y + dTable.yshift;
       dTable.z = dTable.z + dTable.zshift;
    end
end
    

try

switch pars.sortMethod
    case 'byHybe'
        warning('depreciated. Why do you want the data plotted by hyb?');
        % interleave data in the order taken?
        if pars.bins == 0
            % numHybes=max(dTable.hybe);
            numHybes = max(orcaTable.readout(strcmp(orcaTable.dataType,'H')));
        else
           numHybes = pars.bins;
        end
        % initialize the matrices
        polymer = nan(numHybes,length(dims),nPolys);
        distMap = nan(numHybes,numHybes,nPolys);
        spotData = nan(nPolys,2); % just locusX,Y coords
        % populate matrices and compute distances
        for s=1:nPolys
            sTable = dTable(dTable.fs==spotIDs(s),:);
            if ~isempty(sTable)
                    spotData(s,:) = [sTable.locusX(1),sTable.locusY(1)];
                    polymer(sTable.hybe,:,s) = sTable{:,dims}; % Sort by Hyb
                if computeDistMap
                    distMap(:,:,s) = squareform(pdist(polymer(:,:,s)));
                end
            end
        end
        
    case 'byRead'
        % just keep Hybs and Rpts, 
        %       exclude 'B'-barcodes, 'A'-chromatic aligns
        isDat = (strcmp(dTable.dataType,'H') | strcmp(dTable.dataType,'R')) & dTable.readout~=0 ; 
        datTable = dTable(isDat,:);
        nReads = max(datTable.readout(strcmp(datTable.dataType,'H'))); % the readout number of the largest hyb 
        
        % initialize matrices, make space for rpts
        isRpt = strcmp(datTable.dataType,'R');
        rptCodes = unique(datTable.readout(isRpt));
        nRpts = length(rptCodes);
        for r=1:length(rptCodes)
            isR = datTable.readout == rptCodes(r) & isRpt;
            datTable.readout(isR) = nReads + r; % assign new positions to rpt to be at end of matrix
        end
        if pars.bins==0
            bins = nReads+nRpts; 
        elseif pars.bins < nReads + nRpts
            bins = nReads+nRpts; 
        else
            bins = pars.bins; % still need to reassign rpts, even with fixed bins 
        end
        distMap = nan(bins,bins,nPolys);
        polymer = nan(bins,length(dims),nPolys);     
        spotData = nan(nPolys,2);
        
        % populate matrices and compute distances
        for s=1:nPolys    
            sTable = datTable(datTable.fs==spotIDs(s),:);
            if ~isempty(sTable)
                spotData(s,:) = [sTable.locusX(1),sTable.locusY(1)];
                if pars.chromCorrect
                    xyzChromShift = [sTable.xcShift,sTable.ycShift,sTable.zcShift];
                    polymer(sTable.readout,:,s) = sTable{:,dims} + pars.shiftSign*xyzChromShift(:,dims);
                else
                    polymer(sTable.readout,:,s) = sTable{:,dims};
                end 
                if computeDistMap
                    distMap(:,:,s) = squareform(pdist(polymer(:,:,s)));
                end
            end
        end
    otherwise
        error(['sort method ',pars.sortMethod,' not recognized. Valid values: "byHybe" or "byRead"']);
end


catch er
   warning(er.getReport);
   warning('place debug here'); 
end
%%



if pars.dataBasedDriftCorrect
    % compute typical offset
   %  [numSteps,~,nPolys] = size(polymer);
   %  dPoly = nan(numSteps,length(dims),nPolys);  
%     for h=1:numSteps
%        dPoly(h,:,:) = polymer(h,:,:)-polymer(1,:,:); 
%     end
    dPoly = diff(polymer,1);
    meanOff = nanmedian(dPoly(:,:,:),3);
    % new polymer
    newPolymer = polymer;
    newPolymer(2:end,:,:)=newPolymer(2:end,:,:) - repmat(meanOff,[1,1,nPolys]);
    newMap = distMap;
    for n=1:nPolys
        if pars.dims == 3
            newMap(:,:,n) = squareform( pdist( newPolymer(:,:,n)  ));
        else
            newMap(:,:,n) = squareform( pdist( newPolymer(:,:,n)  ));
        end
    end
    % 
     figure(2); clf; 
     subplot(1,3,1); imagesc(nanmedian(distMap,3));
     subplot(1,3,2); imagesc(nanmedian(newMap,3));
     subplot(1,3,3); for d=dims; plot(meanOff(:,d)); hold on; end
    % export
    polymer = newPolymer;
    distMap = newMap;
end


if pars.removeBlank

         % remove blanks
     blnk = isnan(spotData(:,1));
     polymer(:,:,blnk) =[];
     distMap(:,:,blnk) =[];
     spotData(blnk,:) = [];
    
end


