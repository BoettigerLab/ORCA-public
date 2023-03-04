function  [plotMat,coords,max_coords] = ReadWashUMatrix(tableFile,varargin)
% plotMat = ReadWashUMatrix(tableFile) returns a contact matrix from the
% data file specified by tableFile
%
%%  required inputs
% tablefile -- a string which is the complete file path to the saved
% WashUHiC matrix.
% 
%% outputs 
% plotMat - a NxN contact matrix
% coords - the corresponding genome coordinates 
%
% 

defaults = cell(0,3);
defaults(end+1,:) = {'showplot','boolean',false};
defaults(end+1,:) = {'locus','string',''};
defaults(end+1,:) = {'mapRes','nonnegative',0}; % 0 for autodetect
defaults(end+1,:) = {'minReads','nonnegative',0};
defaults(end+1,:) = {'displayRes','nonnegative',0}; % 0 for same as map res
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% convert washU format to Juicebox export format
% wtf washU?
% mixed tabs and commas?  different formats for locus in same file
% .hicup.MAPQ30.KR_k5b.WashU.txt (so it is KR normalized / 'balanced') 
% chrN,start,stop,\t,chrN2:start2-stop2,'\t'
tableIn = readtable(tableFile);
tabData = cellfun(@(x) strsplit(x,'\t'),tableIn{:,3},'UniformOutput',false);
tabData = cat(1,tabData{:});

Var1 = tableIn{:,2}; % this parses fine
[~,Var2,~] = ParseLocusName(tabData(:,2));
Var3 = str2double(tabData(:,3));
newTable = table(Var1,Var2,Var3);

%%
[plotMat,coords,max_coords] = ReadJuiceboxMatrix(newTable,'parameters',pars);
