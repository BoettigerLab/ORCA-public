function WriteDAXFiles(movie, infoFile, varargin)
%--------------------------------------------------------------------------
% WriteDAXFiles(movie, infoFile, varargin)
% This function writes a movie to a .dax file based on the corresponding
% infoFile structure.  
%--------------------------------------------------------------------------
% Inputs:
% movie/int16 array: This array contains the movie data to
%   be saved.  
%
% infoFile/infoFile structure: This array contains the infoFile
% structures that define the properties of each movie.  
%--------------------------------------------------------------------------
% Outputs:
%
%--------------------------------------------------------------------------
% Variable Inputs:
% 
%--------------------------------------------------------------------------
% Alistair Boettiger & Jeffrey Moffitt
% October 3, 2012
% jeffmoffitt@gmail.com
% boettiger.alistair@gmail.com
%
% Version 1.0 - Jeff Moffitt
% Version 1.1 - AB rewritten to prompt for overwrite. 
%    This version is also dependent on the new ParseVariableArguments
%    helper function for building matlab functions.
% Version 1.2
%     Allow a 'saveFullName' to be specified, in place of editing the info
%     file before passing it here.  
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Parse Required Input
%--------------------------------------------------------------------------
if nargin < 2
    error('Both a movie and an infoFile are required');
end

%--------------------------------------------------------------------------
% Parse Variable Input
%--------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'overwrite','boolean',true};
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'saveName','string',''};
defaults(end+1,:) = {'saveFullName','string',''};
defaults(end+1,:) = {'confirmOverwrite','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);
verbose = pars.verbose;


%==========================================================================


if ~isempty(pars.saveFullName)
    [pars.saveFolder,pars.saveName] = fileparts(pars.saveFullName);
end
if ~isempty(pars.saveName)
    pars.saveName = regexprep(pars.saveName,'.dax',''); % saveName should not contain a file end
    infoFile.localName = [pars.saveName,'.dax']; % 
end 
if ~isempty(pars.saveFolder)
    infoFile.localPath = [pars.saveFolder,filesep];
end

% Create dax name
daxName = [infoFile.localName(1:(end-4)) '.dax'];

%--------------------------------------------------------------------------
% Check for if file exists
%--------------------------------------------------------------------------
fileExists = exist([infoFile.localPath daxName],'file');
 
 if fileExists && pars.confirmOverwrite
    disp(['Dax file ',[infoFile.localPath daxName],' exists!']);
    choice = input('overwrite? 1=y,0=n.  ');
    if choice ~= 1 
       pars.overwrite = false; 
    end
 end

if pars.overwrite || ~fileExists
    %----------------------------------------------------------------------
    % Write .ini files
    %----------------------------------------------------------------------
    WriteInfoFiles(infoFile, 'verbose', false);


    %----------------------------------------------------------------------
    % Write DAX file
    %----------------------------------------------------------------------
    if ~isempty( strfind(infoFile.data_type,'little endian') )
        binaryFormat = 'l';
    else
        binaryFormat = 'b';
    end

    fid = fopen([infoFile.localPath daxName], 'w+');
    if fid<0
        error(['Unable to open ' infoFile.localPath daxName]);
    end

    fwrite(fid, ipermute(uint16(movie), [2 1 3]), 'uint16', binaryFormat);

    if verbose
        disp(['Finished writing ' infoFile.localPath daxName]);
    end

    fclose(fid);
end
