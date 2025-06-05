function [imMaxXY,imOutXY,imMaxXZ,imOutXZ] = ShowXY(im,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','handle',100};
defaults(end+1,:) = {'spotXY','cell',[]};
defaults(end+1,:) = {'childHandle','handle',200};
defaults(end+1,:) = {'autoContrastMin','fraction',.3};
defaults(end+1,:) = {'autoContrastMax','fraction',.9999};
defaults(end+1,:) = {'currTime','integer',0};  % current time to show
defaults(end+1,:) = {'showOverlay','boolean',true};
defaults(end+1,:) = {'xylim','array',[]};
defaults(end+1,:) = {'data','struct',[]};
defaults(end+1,:) = {'maxCols','integer',8}; % max number of columns in subplots   
defaults(end+1,:) = {'colormap','colormap',gray(256)};
defaults(end+1,:) = {'colorbar','boolean','false'};
pars = ParseVariableArguments(varargin,defaults,mfilename); 
   


% folderCal2 = 'I:\Jude\2023-10-10_400kb_test\';
% daxFile = [folderCal2,'beads_0001_C1.dax'];
% pars.numChannels =2; 

% 
% pars.currTime  = 1;
% pars.framesPerTime = 10;



%% display XY
numChannels = length(im);

figure(pars.figHandle); clf;
[h,w,z] = size(im{1});

imMaxXY = zeros(h,w,numChannels,class(im{1}));
imOutXY = zeros(h,w,numChannels,class(im{1}));
axs = [];

subplotRows = ceil(numChannels/pars.maxCols);
subplotCols = min(numChannels,pars.maxCols);
syms = {'o','s','>','<','^','v','x','+'};
for c=1:numChannels
    axs(c) = subplot(subplotRows,subplotCols,c);
    imMaxXY(:,:,c) = max(im{c},[],3);
    imOutXY(:,:,c) = IncreaseContrast( imMaxXY(:,:,c),'high',pars.autoContrastMax,'low',pars.autoContrastMin );
    imagesc(imOutXY(:,:,c));
    if ~isempty(pars.spotXY)
        if ~isempty(pars.spotXY{c})
            hold on; plot(pars.spotXY{c}(:,1),pars.spotXY{c}(:,2),syms{c},'MarkerSize',15,'color','r');
        end
    end
    if ~isempty(pars.xylim) % this way we can still zoom out
         xlim(pars.xylim(1:2));  
         ylim(pars.xylim(3:4));
    end
    if pars.colorbar
        colorbar;
    end
end
linkaxes(axs,'xy');
colormap(pars.colormap);

% show overlay as well; 
if numChannels > 1  && pars.showOverlay
    figure(pars.figHandle+1); clf;
    Ncolor(imOutXY);
    if ~isempty(pars.spotXY)
        for c=1:length(pars.spotXY)
            hold on; plot(pars.spotXY{c}(:,1),pars.spotXY{c}(:,2),syms{c},'MarkerSize',15);
        end
    end
    if ~isempty(pars.xylim) % this way we can still zoom out
         xlim(pars.xylim(1:2));  
         ylim(pars.xylim(3:4));
    end
    
    try
      
    % NcolorApp(imOut);  % needs better clean-up routines for open/close
    % behavior
    catch er
        disp(er.getReport);
        disp('debug here')
    end
end

%% display XZ
if false; % z>1
    figure(pars.figHandle+1); clf;
    [h,w,z] = size(im{1});
    imMaxXZ = zeros(z,w,numChannels,class(im{1}));
    imOutXZ = zeros(h,w,numChannels,class(im{1}));
    axsZ = [];
    
        subplotRows = ceil(numChannels/pars.maxCols);
        subplotCols = min(numChannels,pars.maxCols);
    
    for c=1:numChannels
        axsZ(c) = subplot(subplotRows,subplotCols,c);
        imMaxXZ(:,:,c) = max(prermute(im{c},[3,2,1]),[],3);
        imOutXZ(:,:,c) = IncreaseContrast( imMaxXZ(:,:,c),'high',pars.autoContrastMax,'low',pars.autoContrastMin );
        imagesc(imOutXZ(:,:,c));
        if ~isempty(pars.spotXY)
            if ~isempty(pars.spotXY{c})
                try
                hold on; plot(pars.spotXY{c}(:,1),pars.spotXY{c}(:,3),syms{c},'MarkerSize',15,'color','r');
                catch
                end
            end
        end
        if ~isempty(pars.xylim) % this way we can still zoom out
             xlim(pars.xylim(1:2));  
             ylim(pars.xylim(3:4));
        end
        if pars.colorbar
            colorbar;
        end
    end
    linkaxes(axsZ,'xy');
    colormap(pars.colormap);
    
    % show overlay as well; 
    if numChannels > 1  && pars.showOverlay
        figure(pars.figHandle+1); clf;
        Ncolor(imOutXY);
        if ~isempty(pars.spotXY)
            for c=1:length(pars.spotXY)
                hold on; plot(pars.spotXY{c}(:,1),pars.spotXY{c}(:,2),syms{c},'MarkerSize',15);
            end
        end
        if ~isempty(pars.xylim) % this way we can still zoom out
             xlim(pars.xylim(1:2));  
             ylim(pars.xylim(3:4));
        end
        
        try
          
        % NcolorApp(imOut);  % needs better clean-up routines for open/close
        % behavior
        catch er
            disp(er.getReport);
            disp('debug here')
        end
    end
end
