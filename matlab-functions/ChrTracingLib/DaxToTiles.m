function [imTiles,uls,mosaicIms,pars] = DaxToTiles(daxFiles,varargin)
% takes a list of daxfiles (filepaths), loads the dax and the corresponding
% xml files, and creates a mosaic image.  

% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'verbose', 'boolean',true};
defaults(end+1,:) = {'veryverbose', 'boolean',false};
defaults(end+1,:) = {'overwrite', 'boolean',false};
defaults(end+1,:) = {'channel', 'integer',0};
defaults(end+1,:) = {'skipFirst','integer',0};
defaults(end+1,:) = {'figMosaic','integer',20}; % show mosaic
defaults(end+1,:) = {'flipHorizontal','boolean',false};
defaults(end+1,:) = {'flipVertical','boolean',false};
defaults(end+1,:) = {'transpose','boolean',false};
% parameters passed to TilesToMosaic
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'};
defaults(end+1,:) = {'padMosaic','nonnegative',2}; % pad on both sides in multiples of tile size
defaults(end+1,:) = {'downsample','positive',20}; % downsample this fold. 1=no downsampling
defaults(end+1,:) = {'mosaicSize','float',[0,0]}; % if uls have aleady been aligned to some common reference frame, we don't want to disrupt that
% parameters for LoadDax
defaults(end+1,:) = {'driftFolder', 'string',''};  % if empty no drift correct, otherwise, search for an fov_alignTable.csv

pars = ParseVariableArguments(varargin,defaults,mfilename); 
% 

pars.maxProject = true;  
if pars.maxProject
    is3D = false;
else
    is3D = true;
end

nTiles = size(daxFiles,1);
daxProps = LoadDax(daxFiles{1},'justImProps',true,'verbose',pars.veryverbose);
imTiles = cell(nTiles,daxProps.nChns);
uls = zeros(nTiles,2);

for n=1:nTiles
    try
        [imStk,imProps] = LoadDax(daxFiles{n},'parameters',pars,'verbose',pars.veryverbose);
        if isempty(imProps)
            imProps = daxProps;
        end
        for c=1:imProps.chnsLoaded % imProps. = size(chnStk,3);
            if ~is3D
                im = imStk(:,:,c);
            else
                im = imStk(:,:,:,c);
            end
    
    
    %      % camera flips
    %         if imProps.camera.transpose
    %             if ~is3D
    %                 im = im'; % 2D transpose
    %             else
    %                 im = permute(im,[2,1,3]); % 3D transpose
    %             end
    %         end
    %         if imProps.camera.flipVertical
    %             im = flipud(im);
    %         end
    %         if imProps.camera.flipHorizontal
    %             im = fliplr(im);
    %         end
    
    
            % mosaic flips
            %    % orig tr, horz, vert 
            %  21-22, order was tr, vert, horz,  % 4/11/22 -> this doesn't work for scope3.  But T last does. switched  
            %   4/11/22: tr is after horiz.  made tr last. 
    
    
            if imProps.mosaic.flipVertical
                im = flipud(im);
            end   
            if imProps.mosaic.flipHorizontal
                im = fliplr(im);
            end  
    
            if imProps.mosaic.transpose
                if ~is3D
                    im = im'; % 2D transpose
                else
                    im = permute(im,[2,1,3]); % 3D transpose
                end
            end
    
            if pars.verbose
                display(['mosaic flipV ',num2str(imProps.mosaic.flipVertical),...
                        '   mosaic flipH ',num2str(imProps.mosaic.flipHorizontal),...
                        '   mosaic transpose ',num2str(imProps.mosaic.transpose),...
                        '   add. transpose ',num2str(pars.transpose),...
                        '   add. flipV ',num2str(pars.flipVertical),...
                        '   add. flipH ',num2str(pars.flipHorizontal)]);
            end
                    
    
    
            if pars.transpose
                if ~is3D
                    im = im'; % 2D transpose
                else
                    im = permute(im,[2,1,3]); % 3D transpose
                end
            end
            if pars.flipVertical
                im = flipud(im);
            end   
            if pars.flipHorizontal
                im = fliplr(im);
            end   
    
             
    
            imTiles{n,c} = im;
        end
        uls(n,:) = imProps.stageXY(1:2)./imProps.xy2um;

    catch er
        warning(er.getReport);
        disp('debug here');
    end
end


pars.imProps = imProps;

if pars.figMosaic % pars.figMosaic = 20;
    mosaicIms = cell(imProps.chnsLoaded,1);
    for c=1:imProps.chnsLoaded
        mosaicIms{c} = TilesToMosaic(imTiles(:,c),uls,'parameters',pars);
    end
    figure(pars.figMosaic); clf;
    Ncolor( IncreaseContrast(cat(3,mosaicIms{:}),'high',.9999));
    % uls = ulsOut; % update
end