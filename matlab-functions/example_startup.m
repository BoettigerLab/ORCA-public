global matlabFunctionsPath

matlabFunctionsPath = 'C:\Shared\matlab-functions-shared\';
addpath(genpath(matlabFunctionsPath));
% disp('matlabFunctionsPath:');
% disp(matlabFunctionsPath);


global scratchPath
scratchPath = 'F:\Scratch\';
% disp('scratchPath');
% disp(scratchPath);

% this makes it easier to update 
global pyPath
pyPath = 'C:\ProgramData\Anaconda3\';

global condaPrompt;
condaPrompt =  [pyPath,'Scripts\activate.bat && '];

global daoSTORMexe;
daoSTORMexe = [condaPrompt,'cd ', pyPath, 'Lib\site-packages\storm_analysis\ && python.exe ',pyPath,'Lib\site-packages\storm_analysis\daostorm_3d\mufit_analysis.py '];



global defaultXmlFile;  % default fit values for conda
defaultXmlFile = 'C:\Shared\matlab-functions-shared\ChrTracingLib\fitPars_allframe_iters.xml';
% defaultXmlFile = [pyPath,'Lib\site-packages\storm_analysis\test\test_3d_Z.xml'];  % built in test pars   

global cellpose_env;
cellpose_env = [condaPrompt, 'conda activate cellpose2'] ;
% cellpose_env = 'call C:\ProgramData\anaconda3\Scripts\activate.bat  && conda activate cellpose2';
% cellpose_env = 'call C:\ProgramData\anaconda3\Scripts\activate.bat  && conda activate cellpose';

disp('Welcome to Server 2'); 