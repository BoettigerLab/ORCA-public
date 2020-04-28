function ztable = ReadTableFile(filename,varargin)
% read text file tables with arbitrary file endings 
% if it could have been read if it was a txt, it will be read by
% ReadTableFile, which writes a temporary txt copy of the file to the
% scrtach path and tries to read that. 
% filename = 'G:\Alistair\2019-01-05_BeadTests\SlowScan\Hyb_001\ConvZscan_0.off';

global scratchPath
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

defaults = cell(0,3);
defaults(end+1,:) = {'scratchFolder','string',scratchPath};
defaults(end+1,:) = {'overwrite','boolean',true};
defaults(end+1,:) = {'preview','boolean',true};
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[~,object,~] = fileparts(filename);
txtFile = [pars.scratchFolder,object,'.txt'];
if exist(txtFile,'file') && pars.overwrite
    delete(txtFile);
end
[status,cmdout] = system(['copy ',filename,' ',txtFile]); %#ok<ASGLU>  % supresses file copied output
ztable = readtable(txtFile);

if pars.preview
    zmax = min(4,height(ztable));
    disp(ztable(1:zmax,:)); 
end