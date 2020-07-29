function figHandle = SimpleScanReport(localWords, objective, thresholdInds, varargin)
% ------------------------------------------------------------------------
% figHandle = SimpleScanReport(localWords, objective, thresholdInds)
% This function produces a simple visible report on the status of
% ScanThresholds for pipeline 4. 
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
defaults = cell(0,3);
defaults(end+1,:) = {'figHandle', 'handle', []};
defaults(end+1,:) = {'numHybs', 'positive', 16};
defaults(end+1,:) = {'targetOnBits', 'positive', 4};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Create figure handle
% -------------------------------------------------------------------------
try
    figHandle = figure(parameters.figHandle);
catch
    figHandle = figure('Name', 'Scan Report');
end

% -------------------------------------------------------------------------
% Identify current iteraton ind
% -------------------------------------------------------------------------
currentInd = find(sum(thresholdInds,2) == 0, 1, 'first');
if ~isempty(currentInd)
    thresholdInds = thresholdInds(1:(currentInd-1),:);
    objective = objective(1:(currentInd-1));
end

% -------------------------------------------------------------------------
% Plot Objective Function
% -------------------------------------------------------------------------
subplot(2,2,1);
plot(objective);
xlabel('Iteration');
ylabel('Objective');

% -------------------------------------------------------------------------
% Plot Threshold Inds 
% -------------------------------------------------------------------------
subplot(2,2,2);
imagesc(thresholdInds);
xlabel('Hyb');
ylabel('Iteration');

% -------------------------------------------------------------------------
% Compute binary words
% -------------------------------------------------------------------------
binaryWords = de2bi(localWords, parameters.numHybs);
numTargetOnBits = sum(sum(binaryWords,2)==parameters.targetOnBits);

% -------------------------------------------------------------------------
% Plot number of words
% -------------------------------------------------------------------------
subplot(2,2,3);
axesHandle = get(gca);
lines = axesHandle.Children;
if isempty(lines)
    plot(length(localWords)); hold on;
    plot(numTargetOnBits, 'r'); hold on;
else
    data = [numTargetOnBits length(localWords)];
    for i=1:length(lines)
        oldY = get(lines(i), 'YData');
        oldY(end+1) = data(i);
        xData= 1:length(oldY);
        set(lines(i), 'XData', xData, 'YData', oldY);
    end
end
set(gca, 'YScale', 'log');
xlabel('Iteration');
ylabel('Number of Barcodes');

% -------------------------------------------------------------------------
% Plot on bit histogram
% -------------------------------------------------------------------------
subplot(2,2,4);
data = sum(binaryWords,2);
[n, x] = hist(data, 0:max(data));
stairs(x, n);
xlabel('Num On Bits');
ylabel('Counts');

% -------------------------------------------------------------------------
% Force drawing
% -------------------------------------------------------------------------
drawnow;
