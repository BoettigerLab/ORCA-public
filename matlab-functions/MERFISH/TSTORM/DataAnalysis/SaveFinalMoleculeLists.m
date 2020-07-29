function finalSavePaths = SaveFinalMoleculeLists(dataPath, savePath, compositeMListNames, varargin)
% ------------------------------------------------------------------------
% finalSavePaths = SaveFinalMoleculeLists(dataPath, savePath, compositeMListNames, varargin)
% This function writes final molecule lists to savePath from the composite
% lists in dataPath and the thresholdInds.
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% November 28, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
% general parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};

% parameters to read composite bin file name structs
defaults(end+1,:) = {'compositeBinFieldNames', 'cell', {'allCells', 'hybeName', 'hybNum', 'threshold', 'binType'}};
defaults(end+1,:) = {'compositeBinFieldConv', 'cell', {@char, @char, @(x)str2num(x), @str2num, @char}};
defaults(end+1,:) = {'compositeMListType', 'string', 'alist'};

% parameters to write final mList files
defaults(end+1,:) = {'imageTag', 'string', 'STORM'};
defaults(end+1,:) = {'imageMListType', 'string', 'final_alist'};
defaults(end+1,:) = {'notFiducialTag', 'string', 'c1'};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Determine pads for names
% -------------------------------------------------------------------------
numHybs = length(compositeMListNames);
numHybPads = ceil(log10(numHybs));

% -------------------------------------------------------------------------
% Load and save molecule lists
% -------------------------------------------------------------------------
for h=1:length(compositeMListNames)
    % Load molecule list
    mLists = ReadMasterMoleculeList([dataPath compositeMListNames{h}], ...
        'verbose', parameters.verbose, ...
        'compact', true);
    
    % Parse into individual cells and save
    numCells = max([mLists.frame]);
    numCellPads = ceil(log10(double(numCells)));

    for c=1:numCells
        localList = IndexStructure(mLists, [mLists.frame] == c);
        localList.frame = int32(ones(length(localList.frame),1));
        
        fileName = [parameters.imageTag '_' ...
            num2str(h-1, ['%0' num2str(numHybPads) 'd']) '_' ...
            num2str(c-1, ['%0' num2str(numCellPads) 'd']) '_' ...
            parameters.notFiducialTag '_' ...
            parameters.imageMListType ...
            '.bin'];
        try
            WriteMoleculeList(localList, [savePath fileName], 'verbose', parameters.verbose);
        catch
            display('Problem here!');
        end
    end
end
