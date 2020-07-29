function SaveStruct(structIn,varargin)
% SaveStruct(struct);
% SaveStruct(struct, saveFolder);
% SaveStruct(struct,saveFileFullName);

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'requires a structure');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
varInput = varargin(2:end);
parameters = ParseVariableArguments(varInput, defaults, mfilename);


savePath = pwd;
if nargin == 2
    savePath = varargin{1};
end


for f=fieldnames(structIn)';
    eval( [(f{1}), '= {  structIn.',(f{1})  '};'] ) ;
end
structName = inputname(1);

if isdir(savePath)
    saveName = [savePath,filesep,structName,'.mat'];
else
    saveName = savePath; 
end
    
if isempty(structName)
    structName = 'savedStructure';
end
varNames = fieldnames(structIn)';
save(saveName,'structName',varNames{:});

if parameters.verbose
    disp(['wrote ',saveName]);
end

    