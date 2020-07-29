function dataMatrix = CellToMatPad(data,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'padValue','float',nan};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

dataLengths = cellfun(@length,data);
maxObs = max(dataLengths);
nGrps = length(data);
dataMatrix = pars.padValue*ones(maxObs,nGrps);
for n=1:nGrps
    dataMatrix(1:dataLengths(n),n) = data{n};
end
