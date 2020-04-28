 function [spotTable,tforms,fixDist,fsrUsed] = ChromaticAlignFromTable(spotTable,varargin)

% all distance/position values in the input table are now assumed to be 
% in the same units (nanometers). No mixing pixels and nanometers.
% 
% optionally, allows filtering of spots on brightness (ideally this was
% done well the first time). This option also allows explority fits.
% Just specify 'minHeightRef' or 'minHeightOt'.
% 
% ------------------ Notes -------------------------------------
% FOR IMPROVEMENT: we could use some outlier detection -- there is a
% strong prior of how the chromatic shift should look.
% FOR IMPROVEMENT: is polynomial fitting really best? Is there a lower
% dimensional / more accurate neural net representation of chromatic
% distortion?  
% 
% -----------------------------------------------------------------
% Alistair Boettiger
% CC BY Jan 2019


 

defaults = cell(0,3);
% local parameters
defaults(end+1,:) = {'refChn','integer',647};
defaults(end+1,:) = {'minHeightRef','nonnegative',0};
defaults(end+1,:) = {'minHeightOt','nonnegative',0};
% parameters for Polymap3D
defaults(end+1,:) = {'chnColors','colormap',hsv(3)};
defaults(end+1,:) = {'polyOrder','integer',2};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'max2D','nonnegative',350}; % in nm 
defaults(end+1,:) = {'max3D','nonnegative',0}; % in nm
pars = ParseVariableArguments(varargin,defaults,mfilename);


% enforce numerical array for channel
chns = spotTable.chn;
if isstring(chns) || ischar(chns)
    spotTable.chn = cellstr(chns);
    chns = spotTable.chn;
end
if iscell(chns)
    spotTable.chn = cellfun(@str2double,chns);
end

 spotTable.dataType( strcmp(spotTable.dataType,'S')) = 'A';

datChns = unique(spotTable.chn);
otherChns = datChns(datChns ~= pars.refChn);
isAlign = strcmp(spotTable.dataType,'A'); % find align hybes
if sum(isAlign)==0
    error('no hybes designated as "A" or "S" for alignment / color-swap');
end
origReads = unique(spotTable.readout(isAlign)); % find which readouts were imaged in the align hybes
nChns = length(otherChns);

% this is where the chromatic shifts will be saved
spotTable.xcShift = zeros(height(spotTable),1);
spotTable.ycShift = zeros(height(spotTable),1);
spotTable.zcShift = zeros(height(spotTable),1);

tforms = cell(nChns,1);
for c=1:nChns % loop over all non-reference colors
    % Match images of the same spot taken in different colors
    otChn = otherChns(c);
    pOther = cell(length(origReads),1);
    pRef = cell(length(origReads),1);
    for r=1:length(origReads)
        isO = spotTable.readout==origReads(r) & strcmp(spotTable.dataType,'H');
        isA = spotTable.readout==origReads(r) & strcmp(spotTable.dataType,'A');
        isOt = (isA | isO) &  spotTable.chn==otChn; % rpt spots from read r in the 561 chn
        isRef = (isA | isO) & spotTable.chn==pars.refChn; % rpt spots from read r in the 647 chn
        pOther{r} = ...
            [spotTable.x(isOt)+spotTable.locusX(isOt),...
            spotTable.y(isOt)+spotTable.locusY(isOt),...
            spotTable.z(isOt),...
            spotTable.h(isOt),...
            spotTable.fsr(isOt)]; 
            % CantorPair(spotTable.readout(isOt),spotTable.s(isOt))];
        pRef{r} = ...
            [spotTable.x(isRef)+spotTable.locusX(isRef),...
            spotTable.y(isRef)+spotTable.locusY(isRef),...
            spotTable.z(isRef),...
            spotTable.h(isRef),...
            spotTable.fsr(isRef)];
            % CantorPair(spotTable.readout(isRef),spotTable.s(isRef))];
    end
    pOther = cat(1,pOther{:});
    pRef = cat(1,pRef{:});
    % a little filtering 
    pOther(pOther(:,4) < pars.minHeightOt,:) = [];
    pRef(pRef(:,4) < pars.minHeightRef,:) = [];   
    [~,io,ir] = intersect(pOther(:,5),pRef(:,5));
    xyzRef = pRef(ir,1:3);
    xyzWarp = pOther(io,1:3);
    fsrUsed = pRef(ir,5); %fsr indicies
    
    
    % COMPUTE Shift 
    [tform3D,~,fixDist] = Polymap3D(xyzRef,xyzWarp,'parameters',pars);
    tforms{c} = tform3D;
    % add chromatic shifts to data table 
    spotTable = ApplyShiftToTable(spotTable,otChn,tform3D);

    % this is now a function
%     %% ======= Save chromatic shifts to data table =======% 
%     % Add chromatic shifts as xyz shifts for the speicfic points into the table 
%     % (probably best for a new loop over otChns) here, more modular
%     isOt = spotTable.chn == otChn;
%     stageXYZ = [spotTable.x(isOt) + spotTable.locusX(isOt),...
%                 spotTable.y(isOt) + spotTable.locusY(isOt),...
%                 spotTable.z(isOt)];
%     newStageXYZ = tforminv(tform3D,stageXYZ(:,1),stageXYZ(:,2),stageXYZ(:,3));
%     spotTable.xcShift(isOt) = newStageXYZ(:,1) - stageXYZ(:,1);
%     spotTable.ycShift(isOt) = newStageXYZ(:,2) - stageXYZ(:,2);
%     spotTable.zcShift(isOt) = newStageXYZ(:,3) - stageXYZ(:,3);
end 

   