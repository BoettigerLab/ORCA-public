classdef ImageStack < handle
% ------------------------------------------------------------------------
% imDataObj = ImageStack(varargin)
% This class is a storage wrapper for complex stacks of image and meta data. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% August 6, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
% Public properties
properties 
    verbose                 % Determines the verbose state of the class
end

% Visible but protocted properties
properties (SetAccess=protected)
    numFrames               % The number of frames (L)
    frameSize               % The size of a frame (NXM)
    
    indexType               % A list of the meta data types that can be used to index frames
    
    precision               % Record the precision of the provided imageData
    
end

% Hidden properties
properties (Access=private)
    hiddenFieldsToSave = {'index', 'imageData', 'metaData', 'metaDataFormat'}; % A list of the hidden fieds that should be saved/loaded
    
    imageData               % An NxMxL array
    metaData                % A structure of meta data fields with entries for each of the L frames
    
    index                   % A cell array of containers.Map to aid in quick indexing
    
    isMemMapped = false;    % A boolean that determines if the image data is memory mapped
    memMap                  % A memory map to the image data, only used if loaded with memory mapping
    memMapMetaData          % A memory map to the meta data, only used if loaded with memory mapping
    metaDataFormat          % A cell array that describes the organization of the saved meta data
end

% -------------------------------------------------------------------------
% Public Methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = ImageStack(varargin)
        % Create the ImageStack object
        % obj = ImageStack(imageData, metaData)
        % imageData is a NxMxL array
        % metaData is a structure with fields for each type of meta data.
        %    Each field must have L entries of type double or char. 
        
        % -------------------------------------------------------------------------
        % Handle empty class request
        % -------------------------------------------------------------------------
        if nargin < 1
            return;
        else
            imageData = varargin{1};
            metaData = varargin{2};
            varargin = varargin(3:end);
        end

        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3);

        % Parameters for the geometric transformation
        defaults(end+1,:) = {'verbose', 'boolean', false};
        parameters = ParseVariableArguments(varargin, defaults, mfilename);

        % Transfer to object
        fieldsToTransfer = fields(parameters);
        for f=1:length(fieldsToTransfer)
            obj.(fieldsToTransfer{f}) = ...
                parameters.(fieldsToTransfer{f});
        end
        
        % -------------------------------------------------------------------------
        % Save image data
        % -------------------------------------------------------------------------
        obj.imageData = imageData;
        obj.precision = class(imageData);
        
        dataSize = size(imageData);
        
        if length(dataSize) ~= 3
            error('matlabFunctions:invalidArguments', 'imageData must by NxMxL');
        end
        
        % Upate properties of the image stack
        obj.numFrames = dataSize(3);
        obj.frameSize = dataSize(1:2);
        
        % -------------------------------------------------------------------------
        % Check validity of metadata and save
        % -------------------------------------------------------------------------
        if ~isstruct(metaData) & ~strcmp(class(metaData), 'StructureArray')
            error('matlabFunctions:invalidArguments', 'metaData must be a structure');
        end
        obj.metaData = metaData;

        % -------------------------------------------------------------------------
        % Parse meta data to make indices, i.e. look up tables
        % -------------------------------------------------------------------------
        metaDataFields = fields(metaData);
        
        for f=1:length(fields(metaData))
            % Check validity of field
            if length(metaData.(metaDataFields{f})) ~= obj.numFrames
                warning('matlabFunctions:invalidArguments', 'Found invalid meta data. Ignorning');
                continue;
            end

            % Add field name to index types
            obj.indexType{end+1} = metaDataFields{f};
            
            % Find unique entries
            keys = unique(metaData.(metaDataFields{f}));
            
            % Find values (indices of the frames) to match keys
            if length(keys) > 1
                if iscell(keys)
                    values = cellfun(@(x)find(strcmp(metaData.(metaDataFields{f}), x)), keys, 'UniformOutput', false);
                else
                    values = arrayfun(@(x)find(metaData.(metaDataFields{f}) == x), keys, 'UniformOutput', false);
                end
            else
                values = {1:obj.numFrames};
            end
            % Create map
            obj.index{end+1} = containers.Map(keys, values);
        end
        
    end
    
    % -------------------------------------------------------------------------
    % GetImages
    % -------------------------------------------------------------------------
    function [keyValues] = GetKeys(obj, indexType)
        % Return the meta data keys for a given indexType
        % keyValues = obj.GetKeys(indexType)
        
        % -------------------------------------------------------------------------
        % Check input
        % -------------------------------------------------------------------------
        if ~ischar(indexType)
            error('matlabFunctions:invalidInput', 'The index type must be a string');
        end
        
        % -------------------------------------------------------------------------
        % Find index
        % -------------------------------------------------------------------------
        localID = find(strcmp(obj.indexType, indexType));
        
        if isempty(localID)
            error('matlabFunctions:invalidInput', 'Provided index type is invalid.');
        end
        
        % -------------------------------------------------------------------------
        % Return keys
        % -------------------------------------------------------------------------
        localMap = obj.index{localID};
        keyValues = keys(localMap);
      
    end
    % -------------------------------------------------------------------------
    % Export
    % -------------------------------------------------------------------------
    function [imageData, metaData] = Export(obj, exportPath, varargin)
        % Export the image data in a different format. Metadata may or may
        % not be saved.  
        % [imageData, metaData] = obj.GetImages(exportPath, 'format',
        % formatValue);
        
        % -------------------------------------------------------------------------
        % Check required input
        % -------------------------------------------------------------------------
        if nargin < 1 || ~ischar(exportPath)
            error('matlabFunctions:invalidArguments', 'A valid export path must be provided');
        end
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3);

        % Parameters for loading
        defaults(end+1,:) = {'format', {'dax', 'tiff'}, 'dax'};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Switch on format type
        % -------------------------------------------------------------------------
        switch parameters.format
            case 'dax'
                % -------------------------------------------------------------------------
                % Extract file parts
                % -------------------------------------------------------------------------
                [filePath, fileName, fileExt] = fileparts(exportPath);
                
                % -------------------------------------------------------------------------
                % Check provided exportPath against format
                % -------------------------------------------------------------------------
                if ~strcmp(fileExt, '.dax')
                    error('matlabFunctions:invalidArguments', ...
                        'The export path extension does not match the requested format');
                end
                
                % -------------------------------------------------------------------------
                % Prepare information file
                % -------------------------------------------------------------------------
                infFile = CreateInfoFileStructure();
                
                % File name and location
                infFile.localName = [fileName '.inf'];
                infFile.localPath = [filePath filesep()];
                infFile.file = exportPath;
                
                % Image information
                infFile.number_of_frames = obj.numFrames;
                infFile.frame_dimensions = obj.frameSize;
                infFile.frame_size = prod(obj.frameSize);
                infFile.hstart = 1;
                infFile.vstart = 1;
                infFile.hend = obj.frameSize(1);
                infFile.vend = obj.frameSize(2);
                
                % -------------------------------------------------------------------------
                % Write Dax
                % -------------------------------------------------------------------------
                WriteDAXFiles(obj.imageData, infFile, 'verbose', obj.verbose);

            case 'tiff'
                warning('Tiff format is not yet supported.');
        end
            
    end
    % -------------------------------------------------------------------------
    % GetImages
    % -------------------------------------------------------------------------
    function [imageData, metaData] = GetImages(obj, varargin)
        % Slice the image data based on specific values
        % [imageData, metaData] = obj.GetImages('indexType', values);
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        if ~mod(length(varargin)+1,2)
            error('matlabFunctions:invalidArguments', 'Arguments must be passed in pairs');
        end
        
        % -------------------------------------------------------------------------
        % Parse out index tags and values
        % -------------------------------------------------------------------------
        indexTags = varargin(1:2:end);
        keys = varargin(2:2:end);
        
        % -------------------------------------------------------------------------
        % Build boolean index for keeping frames
        % -------------------------------------------------------------------------
        framesToKeep = 1:obj.numFrames;
        
        % -------------------------------------------------------------------------
        % Confirm that all requested index tags exist
        % -------------------------------------------------------------------------
        if ~all(ismember(indexTags, obj.indexType))
            error('matlabFunctions:invalidArguments', ...
                'Some requested index types are not valid');
        end
        
        % -------------------------------------------------------------------------
        % Parse each index tag
        % -------------------------------------------------------------------------
        for i=1:length(indexTags)
            % Find index tag
            localID = strcmp(obj.indexType, indexTags{i});
            
            % Get appropriate map
            localMap = obj.index{localID};
            
            % Make local copy of the desired values
            localKeys = keys{i};
            if ~iscell(localKeys)
                localKeys = {localKeys};
            end
            
            % Get key values
            isKey = localMap.isKey(localKeys);
            
            if ~all(isKey)
                warning('matlabFunctions:invalidKey', ...
                    'Some requested values are not valid');
            end
    
            % Keep only valid keys
            localKeys = localKeys(isKey);
            
            % Compute frames to keep 
            framesToKeepByKey = localMap.values(localKeys);
            newFramesToKeep = [];
            for j=1:length(framesToKeepByKey) % Requests with multiple values should be treated as 'or'
                newFramesToKeep = union(newFramesToKeep, ...
                    framesToKeepByKey{j});
            end
            framesToKeep = intersect(framesToKeep, newFramesToKeep); % Treat as 'and' with different index tags
            
        end
        
        % -------------------------------------------------------------------------
        % Index and return data and meta data
        % -------------------------------------------------------------------------
        if ~obj.isMemMapped % Access already loaded data
            imageData = obj.imageData(:,:, framesToKeep);
        else % image data is memory mapped for efficiency
            imageData = obj.memMap.Data.image(:,:,framesToKeep);
        end

        % Index meta data fields
        if ~obj.isMemMapped
            metaFields = fields(obj.metaData);
            for f=1:length(metaFields)
                metaData.(metaFields{f}) = ...
                    obj.metaData.(metaFields{f})(framesToKeep);
            end
        else
            metaFields = obj.metaDataFormat(:,3);
            for f=1:length(metaFields)
                if strcmp(obj.metaDataFormat{f,1}, 'char') % Char require conversion from uint8
                    metaData.(metaFields{f}) = ...
                        cellstr(char(obj.memMapMetaData.Data.(metaFields{f})(framesToKeep,:)));
                else % Index out of memmap
                    metaData.(metaFields{f}) = ...
                        obj.memMapMetaData.Data.(metaFields{f})(framesToKeep);
                end
            end
        end
    end
    
    % -------------------------------------------------------------------------
    % Overload of length
    % -------------------------------------------------------------------------
    function len = length(obj)
        len = obj.numFrames;
    end
    
    % -------------------------------------------------------------------------
    % Save Function
    % -------------------------------------------------------------------------
    function Save(obj, dirPath)
        % Save the object in a directory specified by dirPath
        % obj.Save(dirPath)
        
        % -------------------------------------------------------------------------
        % Check directory validity
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        [status, ~] = mkdir(dirPath);
        if ~status
            error('matlabFunctions:invalidArguments', 'Invalid directory path');
            return;
        end
        
        % -------------------------------------------------------------------------
        % Define fields to save
        % -------------------------------------------------------------------------
        fieldsToSave = union(properties(obj), obj.hiddenFieldsToSave); % Save both public and private properties
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            PageBreak();
            display(['Saving image stack to ' dirPath]);
            timer = tic;
        end
        
        % -------------------------------------------------------------------------
        % Save fields
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToSave)
            switch fieldsToSave{i}
                case 'imageData' % Save as flat binary file to allow memory mapping
                    if obj.isMemMapped % Simply copy the existing file
                        oldFile = obj.memMap.Filename;
                        copyfile(oldFile, [dirPath 'imageData.bin']);
                    else
                        fid = fopen([dirPath 'imageData.bin'], 'w');
                        fwrite(fid, obj.imageData, obj.precision);
                        fclose(fid);
                    end
                case 'metaData' % Meta date are saved in a format that can be memory mapped
                    if obj.isMemMapped % Simply copy the existing file
                        oldFile = obj.memMapMetaData.Filename;
                        copyfile(oldFile, [dirPath 'metaData.bin']);
                    else
                        % Determine meta data format and organization
                        metaDataFields = fields(obj.metaData);
                        obj.metaDataFormat = cell(length(metaDataFields), 3);
                        obj.metaDataFormat(:,3) = metaDataFields; % Field names

                        % Save fields
                        fid = fopen([dirPath 'metaData.bin'], 'w');
                        for f=1:length(metaDataFields)
                            % Determine type of metaData
                            localData = obj.metaData.(metaDataFields{f});
                            localType = class(localData);
                            obj.metaDataFormat{f, 1} = localType;
                            if strcmp(localType, 'cell') % Only a cell array if the meta data are strings
                                localData = uint8(char(localData)); % Convert to flat array
                                obj.metaDataFormat{f, 1} = 'char'; % Mark as char
                                localType = 'uint8'; % But save as uint8 because memmapfile doesn't handle char
                            end
                            fwrite(fid, localData, localType);

                            % Archive data size
                            obj.metaDataFormat{f, 2} = size(localData);
                        end
                        % Close file pointer
                        fclose(fid);

                    end
                    % Save meta data layout
                    SaveAsByteStream([dirPath 'metaDataFormat.matb'], ...
                        obj.metaDataFormat, 'verbose', obj.verbose);

                otherwise
                    SaveAsByteStream([dirPath fieldsToSave{i} '.matb'], ...
                        obj.(fieldsToSave{i}), 'verbose', obj.verbose);
            end
        end
        
        % -------------------------------------------------------------------------
        % Display progress
        % -------------------------------------------------------------------------
        if obj.verbose
            display(['... completed in ' num2str(toc(timer)) ' s']);
        end
    end
end

% -------------------------------------------------------------------------
% Static methods
% -------------------------------------------------------------------------
methods (Static)
    % -------------------------------------------------------------------------
    % Build a ImageStack object from a saved version
    % -------------------------------------------------------------------------
    function obj = Load(dirPath, varargin)
        % obj = ImageStack.Load(dirPath)
        
        % -------------------------------------------------------------------------
        % Check provided path
        % -------------------------------------------------------------------------
        if dirPath(end) ~= filesep
            dirPath(end+1) = filesep;
        end
        if ~isdir(dirPath)
            error('matlabFunctions:invalidArguments', 'The provided path is not valid');
        end
        
        % -------------------------------------------------------------------------
        % Parse variable inputs
        % -------------------------------------------------------------------------
        defaults = cell(0,3);

        % Parameters for loading
        defaults(end+1,:) = {'verbose', 'boolean', false};
        defaults(end+1,:) = {'memMap', 'boolean', true};
        
        parameters = ParseVariableArguments(varargin, defaults, mfilename);
        
        % -------------------------------------------------------------------------
        % Create empty object (to define fields to load)
        % -------------------------------------------------------------------------
        obj = ImageStack();
        
        % -------------------------------------------------------------------------
        % Define fields to load
        % -------------------------------------------------------------------------
        fieldsToLoad = union(properties(obj), obj.hiddenFieldsToSave);
        fieldsToLoad = setdiff(fieldsToLoad, ...
            {'imageData', 'metaData'}); % These fields are loaded via special means
        
        % -------------------------------------------------------------------------
        % Load properties/data before image
        % -------------------------------------------------------------------------
        for i=1:length(fieldsToLoad)
            switch fieldsToLoad{i}                    
                otherwise
                    obj.(fieldsToLoad{i}) = LoadByteStream([dirPath fieldsToLoad{i} '.matb'], ...
                    'verbose', parameters.verbose);
            end
        end
        
        % -------------------------------------------------------------------------
        % Load image and meta data
        % -------------------------------------------------------------------------
        if ~parameters.memMap
            % Load image data
            fid = fopen([dirPath 'imageData.bin'], 'r');
            obj.imageData = reshape(fread(fid, Inf, [obj.precision  '=>' obj.precision]), [obj.frameSize obj.numFrames]);
            fclose(fid);
            obj.isMemMapped = false;
            
            % Load meta data
            obj.metaDataFormat = LoadByteStream([dirPath 'metaDataFormat.matb'], 'verbose', obj.verbose);
            
            % Open meta data file pointer
            fid = fopen([dirPath 'metaData.bin'], 'r');
            
            % Loop over fields
            metaDataFields = obj.metaDataFormat(:,3);
            for f=1:length(metaDataFields)
                dataSize = obj.metaDataFormat{f,2};
                dataType = obj.metaDataFormat{f,1};
                localData = reshape(...
                    fread(fid, prod(dataSize), [dataType '=>' dataType]), ...
                    dataSize);
                if strcmp(dataType, 'char') % convert char array back to cell array of strings
                    obj.metaData.(metaDataFields{f}) = cellstr(localData);
                else
                    obj.metaData.(metaDataFields{f}) = localData;
                end
            end
            
            % Close file pointer
            fclose(fid);
        else
            % Memory map image data
            obj.memMap = memmapfile([dirPath 'imageData.bin'], ...
                'Writable', false, ...
                'Offset', 0, ...
                'Format', {obj.precision, [obj.frameSize obj.numFrames], 'image'});
            obj.isMemMapped = true;
            
            % Handle char -> uint8 issue: memmapfile does not handle char
            % arrays
            
            charInds = strcmp(obj.metaDataFormat(:,1), 'char');
            localFormat = obj.metaDataFormat;
            localFormat(charInds, 1) = {'uint8'};
            
            obj.memMapMetaData = memmapfile([dirPath 'metaData.bin'], ...
                'Writable', false, ...
                'Offset', 0, ...
                'Format', localFormat);
        end           
    end
    
end % Static methods

end
