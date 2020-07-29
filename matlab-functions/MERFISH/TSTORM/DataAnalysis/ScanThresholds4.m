function [thresholdInds, objective, parameters] = ScanThresholds4(moleculeLists,varargin)
% ------------------------------------------------------------------------
% [bestThetas, report, parameters] = ScanThresholds4(moleculeLists,varargin) 
% This function scans through a series of analyze molecule lists to
% identify the best set of molecules with respect to a provided objective
% function. 
%--------------------------------------------------------------------------
% Necessary Inputs
%--------------------------------------------------------------------------
% Outputs
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger, Jeffrey Moffitt
% November 26, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
% general parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'scanThresholds','boolean',true};
defaults(end+1,:) = {'showFigs','string','on'};
defaults(end+1,:) = {'saveWords', 'boolean', false};
defaults(end+1,:) = {'savePath', 'path', []};

% Parameters for spotfitting
defaults(end+1,:) = {'maxDtoCentroid', 'positive', 1}; 
defaults(end+1,:) = {'minPhotsPerStain', 'positive',0};
defaults(end+1,:) = {'bitOrder','array',[]};

% Parameters for Extract Code And FPKM
defaults(end+1,:) = {'FPKMData', 'struct', []}; % Delimiters for bin files
defaults(end+1,:) = {'codebookPath', 'string', ''}; 
defaults(end+1,:) = {'numHybs', 'positive',[]};
defaults(end+1,:) = {'bitOrder', 'boolean', 1:16}; % Order of bits

% parameters for estimating hyb accuracy
defaults(end+1,:) = {'method', {'minErrorRatePerSpot','minErrorRatePerGene'},'minErrorRatePerGene'};
defaults(end+1,:) = {'oneWeight','nonnegative',0};

% Objective and reporting functions
defaults(end+1,:) = {'objectiveFunction', 'function', []};
defaults(end+1,:) = {'reportFunction', 'function', []};

% Fov selection 
defaults(end+1,:) = {'fovToAnalyze', 'array', []};
defaults(end+1,:) = {'cellsToAnalyze', 'array', []};

defaults(end+1,:) = {'maxIteration', 'positive', 2e4};

% Define initial thresholds
defaults(end+1,:) = {'initalThresholds', 'array', []};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Handle deprecated options
% -------------------------------------------------------------------------
if ~isempty(parameters.cellsToAnalyze)
    parameters.fovToAnalyze = parameters.cellsToAnalyze;
end

% -------------------------------------------------------------------------
% Reorder molecule list to account for imaging order
% -------------------------------------------------------------------------
moleculeLists = moleculeLists(:,parameters.bitOrder,:);
[numCells, numHybes, numThresh] = size(moleculeLists);

if isempty(parameters.fovToAnalyze)
    parameters.fovToAnalyze = 1:numCells;
end

% -------------------------------------------------------------------------
% Initialize thresholds
% -------------------------------------------------------------------------
thresholdInds = zeros(parameters.maxIteration, numHybes);
if isempty(parameters.initalThresholds)
    thresholdInds(1,:) = round(numThresh/2)*ones(1, numHybes);
else
    thresholdInds(1,:) = parameters.initalThresholds;
end

% -------------------------------------------------------------------------
% Iterate over thresholds
% -------------------------------------------------------------------------
iterationNumber = 1;
objective = zeros(1, parameters.maxIteration);
isDone = false;
hybInd = 1;
thresholdInds(1,hybInd) = 1;
while ~isDone
    % -------------------------------------------------------------------------
    % Initialize local words
    % -------------------------------------------------------------------------
    localWords = {};
    for c=1:length(parameters.fovToAnalyze)
        % -------------------------------------------------------------------------
        % Extract MLists
        % -------------------------------------------------------------------------
        mListsByCell = squeeze(moleculeLists(parameters.fovToAnalyze(c),:,:));
        indices = sub2ind(size(mListsByCell), 1:numHybes, thresholdInds(iterationNumber,:));
        mListsByCell = squeeze(mListsByCell(indices));

        % -------------------------------------------------------------------------
        % Calculate Words
        % -------------------------------------------------------------------------
        try
            notEmptyInds = ~cellfun(@isempty, mListsByCell);
            xc = cell(1, length(mListsByCell));
            yc = cell(1, length(mListsByCell));
            xc(notEmptyInds) = cellfun(@(x)x.xc, mListsByCell(notEmptyInds), 'UniformOutput', false);
            yc(notEmptyInds) = cellfun(@(x)x.yc, mListsByCell(notEmptyInds), 'UniformOutput', false);
      
            wordsDetected = FindWords(mListsByCell, [cat(1,xc{:}), cat(1,yc{:})], ...
                'parameters', parameters);
        
            localWords{c} = bi2de(wordsDetected);
        catch
            display('Problems!');
        end
    end
    localWords = cat(1,localWords{:});
    
    % -------------------------------------------------------------------------
    % Save Words
    % -------------------------------------------------------------------------
    if parameters.saveWords & ~isempty(parameters.savePath)
        SaveAsByteStream([savePath 'words_' num2str(iterationNumber) '.matb'], localWords);
    end
    
    % -------------------------------------------------------------------------
    % Evaluate Objective Function
    % -------------------------------------------------------------------------
    objective(iterationNumber) = parameters.objectiveFunction(localWords);
    
    % -------------------------------------------------------------------------
    % Generate Reports
    % -------------------------------------------------------------------------
    if ~isempty(parameters.reportFunction)
        parameters.reportFunction(localWords, objective, thresholdInds);
    end
    
    % -------------------------------------------------------------------------
    % Evaluate Done Condition
    % -------------------------------------------------------------------------
    if iterationNumber >= parameters.maxIteration
        isDone = true;
    end
    
    % -------------------------------------------------------------------------
    % Prepare for next round of iteration
    % -------------------------------------------------------------------------
    if ~isDone
        % Assign next round of indices
        thresholdInds(iterationNumber+1,:) = thresholdInds(iterationNumber,:);

        % Handle the end of a scan case
        if thresholdInds(iterationNumber,hybInd) == numThresh
            scanValues = objective((iterationNumber-numThresh+1):iterationNumber);
            [~, maxInd] = max(scanValues);
            thresholdInds(iterationNumber+1,:) = thresholdInds(iterationNumber,:);
            thresholdInds(iterationNumber+1,hybInd) = maxInd;

            if hybInd == numHybes
                hybInd = 1;
            else
                hybInd = hybInd + 1;
            end
            thresholdInds(iterationNumber+1,hybInd) = 1;
        else
            thresholdInds(iterationNumber+1,hybInd) = thresholdInds(iterationNumber+1,hybInd)+1;
        end

        iterationNumber = iterationNumber + 1;
    end
end

% -------------------------------------------------------------------------
% Save Progress
% -------------------------------------------------------------------------