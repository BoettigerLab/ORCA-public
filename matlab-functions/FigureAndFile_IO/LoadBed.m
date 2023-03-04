function [trk,xs] = LoadBed(bedFile,locusTxt,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'delimiter',{' ','\t'},' '};
defaults(end+1,:) = {'headerLines','integer',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% bedFile = 'U:\GenomeData\Fly\BG3\GSE32749_CTCF_array_BG3_smoothedM.txt';
bedTab = readtable(bedFile,'headerLines',pars.headerLines,'Delimiter',pars.delimiter, 'ReadVariableNames', false);
[chr,st,en] = ParseLocusName(locusTxt);
isChr = strcmp(bedTab{:,1},chr); 
in = bedTab{:,2} > st & bedTab{:,2} < en & isChr;
trk = bedTab{in,4};
xs = bedTab{in,2};