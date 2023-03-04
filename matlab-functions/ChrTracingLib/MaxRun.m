function maxRun = MaxRun(data)

% obs = 6;
% n = round(rand(10,obs))
obs = size(data,2); 
maxRun = zeros(obs,1);
for i=1:obs
measurements = regionprops(logical(data(:,i)), 'Area');
if ~isempty(measurements)
maxRun(i) = max([measurements.Area]);
end
end
