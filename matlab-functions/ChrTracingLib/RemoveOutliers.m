function dataOut = RemoveOutliers(dataIn,varargin)


defaults = cell(0,3); 
defaults(end+1,:) = {'quantile','positive',.95}; % on the scale of the data
defaults(end+1,:) = {'maxBeyond','positive',1.5}; % on the scale of the data
pars = ParseVariableArguments(varargin,defaults,mfilename);

theta = quantile(dataIn,pars.quantile)*pars.maxBeyond;
dataOut = dataIn(dataIn<theta); 