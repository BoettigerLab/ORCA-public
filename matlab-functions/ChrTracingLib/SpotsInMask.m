function fovTable= SpotsInMask(spotXYZH,mask,varargin)
%
%% inputs
% spotXYZH - Nx4xC array of all spots found in the image
%          - can be a table of 
% mask - a bw or labeled array (as from cellpose)
%
% See also 
% ----------------
% SpotsInNucleus
% SpotsInCellsFromDax
% SegmentSpotsPerNucleus


%% Outputs
% fovTable - table of 3DPSF fitted spots
%
%
%  see also SegmentSpotsPerNucleus (which needs only an image)

defaults = cell(0,3);
defaults(end+1,:) = {'maxSpotsPerCell','integer',2};
defaults(end+1,:) = {'troubleshoot','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if istable(spotXYZH) | isstruct(spotXYZH)
    try
        spotXYZH = [Column([spotXYZH.x]),Column([spotXYZH.y]),Column([spotXYZH.z]),Column([spotXYZH.h])];
    catch er
        warning('unrecognized table format. Table must have columns "x","y","z","h"');
        error(er.getReport)
    end
end


spotXY = spotXYZH(:,1:2); 

% assign to nuclei
nCells = size(spotXY,1); 
cellBorders = cell(nCells,1);
fovTable = cell(nCells,1);
for c=1:nCells % c=11
    temp = bwboundaries(mask==c,"noholes"); % sadly this is very slow in a loop
    if ~isempty(temp)
        temp = fliplr(temp{1});  % note the need for column flip
        cellBorders{c} = temp;
        in = inpolygon(spotXY(:,1),spotXY(:,2),...
                       cellBorders{c}(:,1),cellBorders{c}(:,2));
        currSpots = spotXY(in,:);
        if pars.troubleshoot
            figure(4); clf; plot(spotXY(:,1),spotXY(:,2),'r.'); hold on;
            plot(cellBorders{c}(:,1),cellBorders{c}(:,2),'k.-');
        end
        brightness = spotXYZH(in,4);
        z =spotXYZH(in,3);
        numSpots = length(z);
        [~,idx] = sort(brightness,'descend');
        spotsPerCell = min([numSpots,pars.maxSpotsPerCell]);
        selSpots = idx(1:spotsPerCell);
        bright = brightness(selSpots);
        x = currSpots(selSpots,1);
        y = currSpots(selSpots,2);
        z = z(selSpots);
        allele = (1:length(z))';
        
        if pars.troubleshoot
            figure(4);     plot(x,y,'bo'); hold on;
            pause(.01);
        end

        cellNum = ones(spotsPerCell,1)*c;
        cellTable = table(x,y,z,bright,cellNum,allele);
        fovTable{c} = cellTable;
    end
end
fovTable = cat(1,fovTable{:});

