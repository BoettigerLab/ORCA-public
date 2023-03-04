function h = BoxPlotCell(celldata,varargin)
% BoxPlotCell accepts cell array of data of unequal size
% 
%% Inputs
% celldata -- each elment will get its own box in the plot (each is a
% separate set of observations). Cell array may be Nx1 or 1xN
% 
%% optional inputs
% all parameter options for boxplot are allowed, see matlab builtin boxplot


nDatas = length(celldata);
nObs = cellfun(@length,celldata);
maxObs = max(nObs);

matdata = nan(maxObs,nDatas);
for n=1:nDatas
    matdata(1:nObs(n),n) = celldata{n};
end

h = boxplot(matdata,varargin{:}); % need to have the right number of elements 