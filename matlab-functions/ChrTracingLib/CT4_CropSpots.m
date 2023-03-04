function imCrop = CT4_CropSpots(hybIm,imProps,fidTable,varargin)
% crops the whole z-stack in the specified box-size
% note, all box sizes are half width (or radii, not diameter)

defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'align', 'freeType',[]};
defaults(end+1,:) = {'boxSize', 'integer',15};
defaults(end+1,:) = {'cellIDmap', 'freeType',[]}; % cellID map, optional
defaults(end+1,:) = {'chns', 'integer',inf};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% only load requested channels
%    note, we keep an empty slot in the cell array to keep the indexing to
%    macth the channel.
if isinf(pars.chns)
    selectChns = 1:imProps.nChns;
else
    selectChns = pars.chns;
end

chnIm = cell(1,imProps.nChns);
for c=selectChns
    if isempty(pars.align)
        chnIm{c} = hybIm(:,:,:,c);
    else
        chnIm{c} = ApplyReg(hybIm(:,:,:,c),pars.align);  % apply registration
    end
end
% crop by cell
nTraces = height(fidTable);
imCrop = cell(nTraces,imProps.nChns);
for c=selectChns 
    for s=1:nTraces % loop over spots  
        imCrop_z = cell(imProps.zSteps,1);
        if ~isempty(pars.cellIDmap)
            cell_i = fidTable.cellID(s);
        end
        ys = max(1, fidTable.y_pix(s)-pars.boxSize):min(imProps.imSize(2), fidTable.y_pix(s)+pars.boxSize);
        xs = max(1, fidTable.x_pix(s)-pars.boxSize):min(imProps.imSize(1), fidTable.x_pix(s)+pars.boxSize); 
        for z=1:imProps.zSteps                      
            if ~isempty(pars.cellIDmap)
                imMask = pars.cellIDmap(ys,xs);
                im = chnIm{c}(ys,xs,z);  
                im(imMask~=cell_i) = 0;
            else
               im = chnIm{c}(ys,xs,z); 
            end
            imCrop_z{z} = im;
        end
        imCrop{s,c} = cat(3,imCrop_z{:});
    end
end
