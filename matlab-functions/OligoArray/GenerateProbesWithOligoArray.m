function newFilePath = GenerateProbesWithOligoArray(designStruct, varargin)
%--------------------------------------------------------------------------
% newFilePath = GenerateProbesWithOligoArray(varargin)
% This function is an interface between matlab and oligo array. It
% generates an oligo array text file based on the input parameters and
% returns a string to that file path.  
%
%--------------------------------------------------------------------------
% Outputs:
%
% newFilePath/string: A path string to the generated file
%
%--------------------------------------------------------------------------
% Inputs:
%
% probeParameters/struct(required): A structure that specifies the
% parameters of oligoArray
%
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'verbose'/boolean (true): Display or hide function progress
% 'saveDesignStruct'/boolean (true): Save the parameters structure
% 'oligoArrayPath'/string: Path to oligo array executable
% 'blastPath'/string: Path to blast executables
% 'oligoArrayAuxPath'/string: Path to additional oligo array exectutables
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% November 21, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------
requiredFields = {'inputFile', 'databaseFile', 'outputFile', 'failedFile', 'logFile', ...
    'maxOligos', 'minLength', 'maxLength', 'maxDistance', 'minTm', 'maxTm', ...
    'secStructTm', 'crossHybTm', 'minGC', 'maxGC', 'mask', 'minSep', 'numParallel'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global oligoArrayExe;
global legacyBLASTPath;
global oligoArrayAuxPath;
%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
verbose = true;
saveDesignStruct = true;
flags = {'verbose', 'saveDesignStruct', 'oligoArrayExe', 'legacyBLASTPath', ...
    'oligoArrayAuxPath', 'independentProcess'};
independentProcess = false;

%--------------------------------------------------------------------------
% Parse required input
%--------------------------------------------------------------------------
if nargin < 1 || ~isstruct(designStruct)
    error('The first arugment must be a structure');
end

if ~all(ismember(requiredFields, fieldnames(designStruct)))
    errorString = ['The following fields must be included in the experiment structure: '];
    for i=1:length(requiredFields)
        errorString = [errorString requiredFields{i} '\n'];
    end
    error('error:missingFields', errorString);

end
    
%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            case 'saveDesignStruct'
                saveDesignStruct = CheckParameter(parameterValue, 'boolean', 'saveDesignStruct');
            case 'oligoArrayExe'
                oligoArrayExe = CheckParameter(parameterValue,'filePath','oligoArrayExe');
            case 'legacyBLASTPath'
                legacyBLASTPath = CheckParameter(parameterValue, 'fileDir', 'legacyBLASTPath');
            case 'oligoArrayAuxPath'
                oligoArrayAuxPath = CheckParameter(parameterValue, 'fileDir', 'oligoArrayAuxPath');
            case 'independentProcess'
                independentProcess = CheckParameter(parameterValue,'boolean','independentProcess');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Append System Path
%--------------------------------------------------------------------------
setPathCommand = ['PATH=%path%;' legacyBLASTPath ';' oligoArrayAuxPath '& '];
if verbose
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' setPathCommand]);
end

%--------------------------------------------------------------------------
% Generate Command
%--------------------------------------------------------------------------
oligoArrayCommand = ['java -Xmx1024m -jar ' oligoArrayExe ...
    ' -i ' designStruct.inputFile ...                           % Input fasta file
    ' -d ' designStruct.databaseFile ...                        % Database file for comparison
    ' -o ' designStruct.outputFile ...                          % File containing the found probes
    ' -r ' designStruct.failedFile ...                          % File containing the failed probes
    ' -R ' designStruct.logFile ...                             % Log file
    ' -n ' num2str(designStruct.maxOligos) ...                  % The maximum number of oligos per entry
    ' -l ' num2str(designStruct.minLength) ...                  % Minimum oligo length
    ' -L ' num2str(designStruct.maxLength) ...                  % Maximum oligo length
    ' -D ' num2str(designStruct.maxDistance) ...                % Maximum distance between 5' of oligo and 3' of query
    ' -t ' num2str(designStruct.minTm) ...                      % Minimum Tm of oligo
    ' -T ' num2str(designStruct.maxTm) ...                      % Maximum Tm of oligo
    ' -s ' num2str(designStruct.secStructTm) ...                % Temperature for secondary structure prediction
    ' -x ' num2str(designStruct.crossHybTm) ...                 % Threshold for cross-hybridization
    ' -p ' num2str(designStruct.minGC) ...                      % Minimum GC content
    ' -P ' num2str(designStruct.maxGC) ...                      % Maximum GC content
    ' -m ' designStruct.mask ...                                % Sequences features to remove
    ' -g ' num2str(designStruct.minSep) ...                     % Minimum distance between 5' end of oligos
    ' -N ' num2str(designStruct.numParallel)];                  % Number of parallel processors

if independentProcess
    oligoArrayCommand = [oligoArrayCommand ' &'];
end

%--------------------------------------------------------------------------
% Save design structure
%--------------------------------------------------------------------------
if saveDesignStruct
    designStructureFileName = designStruct.outputFile;
    [pathName, fileName] = fileparts(designStructureFileName);
    designStructureFileName = [pathName '\' fileName '_designStruct.mat'];
    if verbose
        display('-----------------------------------------------------------------');
        display('Saving:')
        display(['     ' designStructureFileName]);
    end
    
    save(designStructureFileName, 'designStruct');
end

%--------------------------------------------------------------------------
% Issue Command
%--------------------------------------------------------------------------
if verbose
    display('-----------------------------------------------------------------');
    display('Issuing:')
    display(['     ' oligoArrayCommand]);
end
system([setPathCommand oligoArrayCommand]);

