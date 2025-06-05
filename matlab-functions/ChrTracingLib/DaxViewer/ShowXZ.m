function [imMaxZ,imOutZ] = ShowXZ(im,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'childHandle','handle',200};
defaults(end+1,:) = {'autoContrastMin','fraction',.3};
defaults(end+1,:) = {'autoContrastMax','fraction',.9999};
defaults(end+1,:) = {'showOverlay','boolean',true};
defaults(end+1,:) = {'xylim','array',[]};
defaults(end+1,:) = {'colormap','colormap','gray'};
defaults(end+1,:) = {'colorbar','boolean','true'};
pars = ParseVariableArguments(varargin,defaults,mfilename); 
   


%% show the xz projection for the current view
figure(pars.childHandle); clf;
numChannels = length(im);

[h,w,z] = size(im{1});
if ~isempty(pars.xylim)
    xl = round(pars.xylim(1:2)); 
    yl = round(pars.xylim(3:4));
    w = xl(2)-xl(1)+1;
else % use entire image (typically cropped before passing)
    % xl = get(gca,'xlim');
    % yl = get(gca,'ylim');
    xl = round([1,w]);
    yl = round([1,h]);
end   
    imMaxZ = zeros(z,w,numChannels,class(im{1}));
    imOutZ = zeros(z,w,numChannels,class(im{1}));
    axs = [];
    for c=1:numChannels
        axs(c) = subplot(1,numChannels,c);
        try
        imMaxZ(:,:,c) = squeeze(max(  permute(im{c}(yl(1):yl(2),xl(1):xl(2),:),[3,2,1])     ,[],3));
        imOutZ(:,:,c) = IncreaseContrast( imMaxZ(:,:,c) ); % 'high',pars.autoContrastMax,'low',pars.autoContrastMin
        imagesc(imOutZ(:,:,c));

        if pars.colorbar
            colorbar;
        end
        catch er
            warning(er.getReport)
            disp('debug')
        end
    end
    linkaxes(axs,'xy');
    colormap(pars.colormap);



% show overlay as well; 
if numChannels > 1
    figure(pars.childHandle+1); clf;
    Ncolor(imOutZ);
end
