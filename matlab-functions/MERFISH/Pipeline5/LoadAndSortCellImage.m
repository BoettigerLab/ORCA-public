function alignedIm = LoadAndSortCellImage(cellName,idxBits,c,varargin)
% Child function of Decode Pixels (step 2)
global figureSavePath;

defaults = cell(0,3);
defaults(end+1,:) = {'borderSize', 'nonnegative', 45}; %  % down scale this if we drop upsampling.
defaults(end+1,:) = {'minPixels', 'nonnegative', 8}; %  % down scale this if we drop upsampling. 
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'moreVerbose', 'boolean', true}; % 
defaults(end+1,:) = {'overwrite', 'boolean', true}; % 
defaults(end+1,:) = {'savePath', 'string', figureSavePath}; % 

parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main function
alignedIm = ReadDax(cellName,'verbose',false); % readDax defaults to double 
alignedIm = alignedIm(:,:,idxBits); 
numBits = size(alignedIm,3);

% remove the border
alignedIm(1:parameters.borderSize,:,:) = [];
alignedIm(end-parameters.borderSize+1:end,:,:) = []; % important
alignedIm(:,1:parameters.borderSize,:) = []; % important
alignedIm(:,end-parameters.borderSize+1:end,:) = [];


% just plotting
if parameters.saveFigures || strcmp(parameters.figVis,'on');
    normIm = figure('Name','normIm','visible',parameters.figVis);  
    Ncolor(uint16(alignedIm),hsv(numBits)); FluorImage('trim',0); 
    set(gcf,'Units','Inches','Position',[4 1 8 8]);
    SaveFigure(normIm,'name',['normIm_',num2str(c,'%02d')],'formats',{'png'},...
        'overwrite',parameters.overwrite,'saveData',parameters.saveFigures,...
        'closeFig',parameters.closeOnSave,'verbose',parameters.moreVerbose,...
        'savePath',parameters.savePath);
    colordef white;
end
