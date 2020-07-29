function [overlays,validation] = ChrTracer3p2_ValidateDriftFix(fidFrames,varargin)

% Global variables
global figureSavePath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'eTable','freeType',[]};
defaults(end+1,:) = {'selectFOVs','freeType',inf}; %

% summary plot specific parameters
defaults(end+1,:) = {'contrastLow','fraction',.5};
defaults(end+1,:) = {'contrastHigh','fraction',.9995};
defaults(end+1,:) = {'recordValidation','boolean',true};
defaults(end+1,:) = {'saveFigFile', 'boolean', true}; 
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'saveFormats', 'cell',{'fig','png'}}; %  {'fig','png'}, {}

% % obsolete
% defaults(end+1,:) = {'saveData', 'boolean', false}; 
% defaults(end+1,:) = {'gain', 'positive', 1}; 
% defaults(end+1,:) = {'showFidXY','boolean',true};
% defaults(end+1,:) = {'rescaleTile','positive',.1};
% defaults(end+1,:) = {'colormap','colormap','jet'};


pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main function

if isempty(pars.saveFormats)
    pars.saveData = false;
end

numFOVs = length(fidFrames); 
selectFOVs = true(1,numFOVs);
if any(pars.selectFOVs) && ~any(isinf(pars.selectFOVs))
    selectFOVs  = ~selectFOVs;
    selectFOVs(pars.selectFOVs) = true;
end
fovToRun = selectFOVs;


% first compute all the overlays, which takes some time, then plot them
if pars.verbose
    disp('assemblying overlays, please wait...');
end
validation = true(numFOVs,1);
overlays = cell(numFOVs,1);
for f = find(selectFOVs)
    fiducialAlignFrames = fidFrames{f};
    [overlays{f},fovToRun(f)] = SaveFidAlignFrames(fiducialAlignFrames,'fov',f,'parameters',pars);
end

if pars.verbose
    disp('displaying overlays for validation...');
end

figHandle = [];
validation(~fovToRun) = false; % don't prompt repeating things we have no info on.
for f=find(fovToRun)
    if ~isempty(figHandle)
        delete(figHandle);
    end
    
    % xy overlay
    figName = ['fov',num2str(f,'%03d'),'_fid_overlayFig'];
    if exist([figureSavePath,figName,'.fig'],'file')
        uiopen([figureSavePath,figName,'.fig'],1);
        figHandle = gcf;
    else
        figHandle = figure(1); clf; 
        imagesc(overlays{f});
        figHandle.Name = figName;
    end
    
    % zoom figure
    % [to add later]
        
    if pars.recordValidation
        c1 = uicontrol(figHandle,'Style','radiobutton','String','Bad?');
        uicontrol(figHandle,'Position',[20,00,60,20],'Style','pushbutton','String','Next','callback',@CloseWindow);
        uiwait;
        validation(f) = c1.Value;
        delete(figHandle);
    end
   pause(.01); 
end

redo = find(validation);
if ~isempty(redo)
    disp('Return to DriftCorrection and use "step pars" to select the following FOVs for realignment: ')
    disp(num2str(redo));
else
    disp('All alignments look good! proceed to Next Step');  
end



function CloseWindow(src,event)
uiresume;
