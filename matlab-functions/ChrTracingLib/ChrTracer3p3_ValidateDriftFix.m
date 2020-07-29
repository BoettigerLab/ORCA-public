function validation = ChrTracer3p3_ValidateDriftFix(saveFolder,varargin)

% Global variables
global figureSavePath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'selectFOVs','freeType',inf}; %
defaults(end+1,:) = {'recordValidation','boolean',true};

pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main function

overlayNames = ls([saveFolder,'*_fid_overlayFig.fig']);
overlayFiles = cellstr(overlayNames);
if ~isempty(overlayFiles)
    overlayFiles = strcat(saveFolder,overlayFiles); %#ok<NASGU>
    fovWithData = str2num(overlayNames(:,4:6))'; %#ok<ST2NM>
else
   fovWithData = [];
end

if ~isinf(pars.selectFOVs)
    fovToRun = intersect(pars.selectFOVs,fovWithData);
else
    fovToRun = fovWithData;
end

numFOVs = length(fovToRun); 
validation = false(numFOVs,1);
figHandle = [];
for f=fovToRun % fovToRun is not boolean
    if ~isempty(figHandle)
        delete(figHandle);
    end
    
    % xy overlay
    figName = ['fov',num2str(f,'%03d'),'_fid_overlayFig'];
    if exist([figureSavePath,figName,'.fig'],'file')
        uiopen([figureSavePath,figName,'.fig'],1);
        figHandle = gcf;
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

redo = fovToRun(validation);
if ~isempty(redo)
    disp('Returning to Global Drift Correction. Adjust "Step Pars" to fix the alignment of the following FOVs: ')
    disp(num2str(redo));
else
    disp('All alignments look good! proceed to Next Step');  
end



function CloseWindow(src,event)
uiresume;
