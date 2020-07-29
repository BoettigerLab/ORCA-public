function [thresholdInds, objective, parameters] = ScanThresholds4p1(moleculeLists,varargin)
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
defaults(end+1,:) = {'numHybs', 'positive',[]};
defaults(end+1,:) = {'bitOrder', 'boolean', 1:16}; % Order of bits

% Objective and reporting functions
defaults(end+1,:) = {'objectiveFunction', 'function', @Default};
defaults(end+1,:) = {'reportFunction', 'function', []};
defaults(end+1,:) = {'reportFunction2', 'function', []};
defaults(end+1,:) = {'reportFunction3', 'function', []};
defaults(end+1,:) = {'reportFunction4', 'function', []};

% Cell Inds
defaults(end+1,:) = {'cellsToAnalyze', 'array', []};
defaults(end+1,:) = {'maxIteration', 'positive', 1e3};
defaults(end+1,:) = {'startingThresholdInd', 'positive', 5};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Reorder molecule list to account for imaging order
% -------------------------------------------------------------------------
moleculeLists = moleculeLists(:,parameters.bitOrder,:);
[numCells, numHybes, numThresh] = size(moleculeLists);

if isempty(parameters.cellsToAnalyze)
    parameters.cellsToAnalyze = 1:numCells;
end

% -------------------------------------------------------------------------
% Iterate over thresholds
% -------------------------------------------------------------------------
iterationNumber = 1;
thresholdInds = zeros(parameters.maxIteration, numHybes);
objective = zeros(1, parameters.maxIteration);
maxAboveBlank = zeros(1, parameters.maxIteration);
ratio4to1 = zeros(1, parameters.maxIteration);
fpkmCorr = zeros(1, parameters.maxIteration);
maxHSRatio  = zeros(1, parameters.maxIteration);

thresholdInds(1,:) = parameters.startingThresholdInd*ones(1, numHybes);
isDone = false;
hybInd = 1;
thresholdInds(1,hybInd) = 1;
while ~isDone
    % -------------------------------------------------------------------------
    % Initialize local words
    % -------------------------------------------------------------------------
    localWords = {};
    for c=1:length(parameters.cellsToAnalyze)
        % -------------------------------------------------------------------------
        % Extract MLists
        % -------------------------------------------------------------------------
        mListsByCell = squeeze(moleculeLists(parameters.cellsToAnalyze(c),:,:));
        indices = sub2ind(size(mListsByCell), 1:numHybes, thresholdInds(iterationNumber,:));
        mListsByCell = squeeze(mListsByCell(indices));

        % -------------------------------------------------------------------------
        % Calculate Words
        % -------------------------------------------------------------------------
        mLists = [mListsByCell{:}];
        wordsDetected = FindWords(mListsByCell, [cat(1,mLists.xc), cat(1,mLists.yc)], ...
            'parameters', parameters);
        
        localWords{c} = bi2de(wordsDetected);
        
    end
    localWords = cat(1,localWords{:});
    
    % -------------------------------------------------------------------------
    % Save Words
    % -------------------------------------------------------------------------
    if parameters.saveWords & ~isempty(parameters.savePath)
        SaveAsByteStream([parameters.savePath 'words_' num2str(iterationNumber) '.matb'], localWords);
    end
    % -------------------------------------------------------------------------
    % Evaluate Objective Function
    % -------------------------------------------------------------------------
    objective(iterationNumber) = parameters.objectiveFunction(localWords);
    
    
    % -------------------------------------------------------------------------
    % Determine New Thresholds
    % -------------------------------------------------------------------------
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
    % -------------------------------------------------------------------------
    % Generate Reports
    % -------------------------------------------------------------------------
    if ~isempty(parameters.reportFunction)
        maxAboveBlank(iterationNumber) = parameters.reportFunction(localWords);
        fpkmCorr(iterationNumber) = parameters.reportFunction2(localWords); 
        ratio4to1(iterationNumber) = parameters.reportFunction3(localWords); 
        maxHSRatio(iterationNumber)= parameters.reportFunction4(localWords); 
    end
    
    % -------------------------------------------------------------------------
    % Evaluate Done Condition
    % -------------------------------------------------------------------------
    if iterationNumber >= parameters.maxIteration
        isDone = true;
    end
    if ~mod(iterationNumber,20)
        display(iterationNumber);
        maxHSRatio((iterationNumber-19):iterationNumber)
        fpkmCorr((iterationNumber-19):iterationNumber)
    end
    iterationNumber = iterationNumber + 1;
end

% -------------------------------------------------------------------------
% Save Progress
% -------------------------------------------------------------------------

parameters.maxAboveBlank = maxAboveBlank;
parameters.fpkmCorr = fpkmCorr; 
parameters.ratio4to1 = ratio4to1;
parameters.maxHSRatio = maxHSRatio;


