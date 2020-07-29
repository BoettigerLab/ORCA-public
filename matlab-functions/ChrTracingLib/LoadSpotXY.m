function regData = LoadSpotXY(saveFolder,varargin)
% find regData in saveFolder for fov requested (defaults to 1)
% name format is fixed/protected: fovNNN_regData.csv

defaults = cell(0,3);
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'suffix','string','_selectSpots.csv'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

tableName = [saveFolder,'fov',num2str(pars.fov,'%03d'),pars.suffix];
if exist(tableName,'file') 
       tableData = readtable(tableName);
       regData = tableData{:,:};
else
    regData = [];
end