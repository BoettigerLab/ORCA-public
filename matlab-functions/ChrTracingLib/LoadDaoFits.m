function data = LoadDaoFits(binName,varargin)
% load either a .hdf5 or a .csv file created by storm_analysis (e.g.
% DaoSTORM).  

global condaPrompt pyPath  %#ok<GVMIS> % load python paths from startup; 

% Parse default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.verbose
    disp(['loading ',binName,' ... ']);
end

tic
csvFile = regexprep(binName,'.hdf5','.csv');
% if a csv-version of the data has not already been created, generate
% it using the storm_analysis converter:
if ~exist(csvFile,'file')
    pythonSA = [condaPrompt,'cd ', pyPath, 'Lib\site-packages\storm_analysis\ && python.exe '];
    hd5_to_txt = './sa_utilities/hdf5_to_txt.py';
    cmdOut = [pythonSA,hd5_to_txt ' --hdf5 ',binName,' --txt ',csvFile];
    if pars.verbose
        disp(cmdOut);
    end
    system(cmdOut);
end
% now we can load the table.
data = readtable(csvFile); 
t = toc;
if pars.verbose
    disp(['data loaded in ',num2str(t,4),' seconds']);
end
   