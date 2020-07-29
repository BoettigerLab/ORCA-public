function estimates = SeparationPlane(setA, setB)
% find the plane that is maxidistant to separate the data points
% setA and setB are m x n matrices where m is the number of samples and n
% is the number of features measured for each sample
% separationPlane([1 5; 2 3; 5 1;], [6 2; 3 5; 7 7; 8 3;]);

% output plane is of the form ax + by + cz + ... + d = 0

% distance from a point (x, y, z) is:
% |ax + by + cz + d| / sqrt(a^2 + b^2 + c^2)
% bottom part obviously disappears if a,b,c normalized to each other

if nargin < 2 || size(setA, 2) ~= size(setB, 2) || ~isnumeric(setA) || ~isnumeric(setB)
    error('Input must be two numeric matrices with the same second dimension')
end

% define a display variable
showDisplay = ~nargout;

% normalize the dimensions so that they contribute equally
normCoeff = [min([setA; setB], [], 1); max([setA; setB], [], 1) - min([setA; setB], [], 1)]';
normData = ([setA; setB] - repmat(normCoeff(:, 1)', size([setA; setB], 1), 1)) ./ repmat(normCoeff(:, 2)', size([setA; setB], 1), 1);
setA = normData(1:size(setA, 1), :) * 10;
setB = normData(size(setA, 1) + 1:end, :) * 10;

% show a projection plane during fitting if no output is requested
if showDisplay
    figure('numbertitle', 'off'), plot(setA(:, 1), setA(:, 2), 'lineStyle', 'none', 'marker', '+', 'markeredge', [1 0 0]);
    hold on
    plot(setB(:, 1), setB(:, 2), 'lineStyle', 'none', 'marker', '+', 'markeredge', [0 0 1]);
    lineHandle = line(1,1, 'color', [0 0 0]);
    set(gca, 'xlim', [-.5 10.5], 'ylim', [-.5 10.5]);
end

centroidA = mean(setA);
centroidB = mean(setB);

if showDisplay
    plot(centroidA(1), centroidA(2), 'marker', 'o', 'color', 'r')
    plot(centroidB(1), centroidB(2), 'marker', 'o')
end

setA(:, end + 1) = 1;
setB(:, end + 1) = 1;
% Call fminsearch with the centroid perpendicular bissector as a seed
estimates = fminsearch(@minFun, [(centroidB - centroidA) / sqrt(sum((centroidB - centroidA).^2)) -sqrt(sum(mean([centroidA; centroidB]).^2))], optimset('MaxFunEvals', 1000000, 'Display', 'none'));

% denormalize the estimates
normEstimates = estimates;
estimates(1:end - 1) = estimates(1:end - 1) ./ (normCoeff(:, 2)' ./ 10);
estimates(end) = -estimates(1:end - 1) * ((-normEstimates(end) ./ sum(normEstimates(1:end - 1).^2) * normEstimates(1:end - 1) .* (normCoeff(:,2)' ./ 10) + normCoeff(:,1)')');

% make plane vector unit
estimates = estimates ./ sqrt(sum(estimates(1:end - 1).^2));

if showDisplay
    outString = '';
    for i = 1:numel(estimates)
        outString = [outString ', ' char(64 + i) ' = ' sprintf('%1.3f', estimates(i))];
    end
    set(gcf, 'name', [get(gcf, 'name') outString]);
end

    function errorSum = minFun(params)
        normFactor = sqrt(sum(params(1:end - 1).^2)); % normalization factor
        % sum of the distance of each point from the plane, with positive distance in the direction of the set's centroid
        errorSum = sum(distanceFunction(([sign(params(1:end - 1) * (centroidB - centroidA)') .* params(1:end - 1) params(end)] * setA') ./ normFactor)) + sum(distanceFunction(([sign(params(1:end - 1) * (centroidB - centroidA)') .* params(1:end - 1) params(end)] * setB') ./ -normFactor));
        
        % show fit line
        if showDisplay
            set(lineHandle, 'xdata', [min([setA(:, 1); setB(:, 1)]) max([setA(:, 1); setB(:, 1)])], 'ydata', (-params(3) - params(1) .* [min([setA(:, 1); setB(:, 1)]) max([setA(:, 1); setB(:, 1)])]) ./ params(2));
            set(gcf, 'name', ['Error = ' sprintf('%1.1f', errorSum)]);
            pause(.1);
        end
        
        function distance = distanceFunction(inData)
            % an exponential function of distance with a step of 10000 when
            % crossing from positive to negative distance
            distance = exp(inData) + (inData > 0) .* 10000;
        end        
    end
end