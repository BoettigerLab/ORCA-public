function [filePaths,fileNames] = FindFileInSubFolders(parentFolder,fileName,varargin)
% ----------------------------------------------------------------------- % 
% returns a list of complete file paths to all files matching "fileName"
% (wildcard * characters are allowed). 
% 
% Example 1: filePaths = FindFileInSubFolders('C:\Data\','*.i4d')
% returns all .i4d files in 'C:\Data\' or any of its subfolders
% 
% Example 2: Return just the file names, not the full path
% [~,fileNames] = FindFileInSubFolders('C:\Data\','*.i4d')
% 
% [~,foundSpotTables] = FindFileInSubFolders(saveFolder,'fov*AllFits.csv')
% ----------------------------------------------------------------------- % 
% Alistair Boettiger
% CC BY NC Aug 9 2017

% fileName = '*AllFits.csv'
% fileName = '*.csv'

defaults = cell(0,3);
defaults(end+1,:) = {'depth','nonnegative',inf};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.depth == 0
    fileNames = cellstr(ls([parentFolder,filesep,fileName]));
    filePaths = strcat(parentFolder,fileNames);
else

    allFolders = strsplit(genpath(parentFolder),';');
    allFolders = allFolders(1:end-1); % remove terminal ';'
    nFolders = length(allFolders);
    fullPaths = cell(nFolders,1); 
    fileNames = cell(nFolders,1); 
    for i=1:nFolders
        fileNames{i} = cellstr(ls([allFolders{i},filesep,fileName]));
        fullPaths{i} = strcat([allFolders{i},filesep], fileNames{i});
    end
    filePaths = cat(1,fullPaths{:});
    fileNames = cat(1,fileNames{:});
    filePaths(cellfun(@isempty,filePaths)) = [];
    fileNames(cellfun(@isempty,fileNames)) = [];
    [~,~,fileType] = fileparts([parentFolder,fileName]);
    keep = contains(filePaths,fileType);
    filePaths = filePaths(keep);


end
