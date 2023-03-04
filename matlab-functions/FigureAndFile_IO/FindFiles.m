function foundFiles = FindFiles(strIn,varargin)
% a shortcut for cellstr( ls( strIn ));
% which is smart enough to return an empty cell rather than a 1x1 cell of
%     an empty character array, if nothing is found, rather than an empty 
%     cell.
% also allows collection of folders only, enforcing the convention that
% folders end with a filesep. 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'fullPath','boolean',true};
defaults(end+1,:) = {'removeDot','boolean',true}; % file up dir "." and".."    
defaults(end+1,:) = {'onlyFolders','boolean',false};
defaults(end+1,:) = {'onlyFiles','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[folder,fname,ftype] = fileparts(strIn);
if pars.onlyFolders    
    dirOut = dir(strIn);
    if pars.fullPath
        foundFiles = strcat([folder,filesep],{dirOut([dirOut.isdir]).name},filesep)';
    else
        foundFiles = strcat({dirOut([dirOut.isdir]).name},filesep)';
    end
else
    
    foundFiles = ls(strIn);
    if ~isempty(foundFiles) 
        if pars.fullPath
            foundFiles = strcat([folder,filesep],cellstr(foundFiles));
        else
            foundFiles = cellstr(foundFiles);
        end
    else
        foundFiles = {};
    end 
end

if pars.removeDot
    if pars.fullPath
        toRemove = contains(foundFiles,'\.');
        foundFiles(toRemove) = [];
    else
        toRemove = strcmp(foundFiles,'.\');
        toRemove = toRemove | strcmp(foundFiles,'..\');
        toRemove = toRemove | strcmp(foundFiles,'.');
        toRemove = toRemove | strcmp(foundFiles,'..');
        foundFiles(toRemove) = [];
    end
end
