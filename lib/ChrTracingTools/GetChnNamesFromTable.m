function [fidChn,dataChns] = GetChnNamesFromTable(eTable,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'hyb', 'integer', 1};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin, defaults, mfilename);


h = pars.hyb;
channels = strsplit(regexprep(eTable.channels{h},'[^A-Za-z0-9,]',''),',');
nBufFrames = eTable.bufferFrames(h);
totFrames  = eTable.totalFrames(h);
frameChannels = cell(totFrames,1);
numChns = length(channels);
fidChannel = strcmp(channels,num2str(eTable.fiducialChannel(h)));
dataChns = channels(~fidChannel);
fidChn = channels(fidChannel);