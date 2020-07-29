classdef OTMap2 < handle
% ------------------------------------------------------------------------
% otMap = OTMap2(initialData, varargin)
% This class provides an interface to a very fast key/value interface.
%--------------------------------------------------------------------------
% Necessary Inputs
% initialData -- An 2xN array of key (1) and value (2) pairs. Both must be
% doubles. The key values do not need to be unique
%--------------------------------------------------------------------------
% Methods
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% None
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 20, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define properties
% -------------------------------------------------------------------------
properties (GetAccess=private)
    data
end

% -------------------------------------------------------------------------
% Define methods
% -------------------------------------------------------------------------
methods
    % -------------------------------------------------------------------------
    % Define constructor
    % -------------------------------------------------------------------------
    function obj = OTMap2(initialData)
        % -------------------------------------------------------------------------
        % Check intput
        % -------------------------------------------------------------------------
        if nargin < 1 || isempty(initialData)
            obj.data = containers.Map('KeyType','double','ValueType', 'double');
            return
        end
        if ~isa(initialData, 'double') || size(initialData,1) ~= 2
            error('matlabFunctions:invalidArguments', ...
                ['initialData must be a double of size 2xN']);
        end
        
        % -------------------------------------------------------------------------
        % Find unique keys and accumulate values
        % -------------------------------------------------------------------------         
        [uniqueKeys, ~, ic] = unique(initialData(1,:));
        values = accumarray(ic, initialData(2,:), [])';
        
        % -------------------------------------------------------------------------
        % Set data
        % -------------------------------------------------------------------------         
        obj.data = containers.Map(uniqueKeys, values);
    end
    
    % -------------------------------------------------------------------------
    % AddToMap
    % -------------------------------------------------------------------------
    function AddToMap(obj, newData)
        % -------------------------------------------------------------------------
        % Check data
        % -------------------------------------------------------------------------         
        if ~isa(newData, 'double') || size(newData,1) ~=2
            error('Invalid data format');
        end
        
        % -------------------------------------------------------------------------
        % Find unique entries and sum as needed
        % -------------------------------------------------------------------------
        [keys, ~, ic] = unique(newData(1,:));
        values = accumarray(ic, newData(2,:), [])';
        
        % -------------------------------------------------------------------------
        % Find overlapping values, sum where needed, and reassign
        % -------------------------------------------------------------------------
        validKeys = obj.data.isKey(num2cell(keys));
        if any(validKeys)
            updatedValues = obj.data.values(num2cell(keys(validKeys)));
            updatedValues = [updatedValues{:}] + values(validKeys);
            oldKeys = keys(validKeys);
        else
            oldKeys = [];
            updatedValues = [];
        end
        newKeys = keys(~validKeys);
        newValues = values(~validKeys);
        newMap = containers.Map([oldKeys newKeys], [updatedValues newValues]);
        
        obj.data = [obj.data; newMap];
    end
    
    % -------------------------------------------------------------------------
    % SubtractFromMap
    % -------------------------------------------------------------------------
%     function obj = SubtractFromMap(obj, newData)
%         warning('Not yet implemented!')
%         % -------------------------------------------------------------------------
%         % Check data
%         % -------------------------------------------------------------------------         
%         if ~isa(newData, 'double') || size(newData,1) ~=2
%             error('Invalid data format');
%         end
%         
%         % -------------------------------------------------------------------------
%         % Identify keys to subtract: issue warning if keys are not in both
%         % maps
%         % -------------------------------------------------------------------------         
%         validKeys = obj.data.isKey(newData(1,:));
%         if ~all(validKeys)
%             warning('matlabFunctions:missingKeys', ...
%                 'Some keys are present for subtraction and will be ignored!');
%         end
% 
%         % -------------------------------------------------------------------------
%         % Subtract values and reassign to the keys
%         % -------------------------------------------------------------------------         
%         if any(validKeys)
%             updatedValues = obj.data.values(num2cell(newData(1,validKeys)));
%             updatedValues = [updatedValues{:}] - newData(2,validKeys);
%             oldKeys = keys(validKeys);
%             
%             %% UNDER CONSTRUCTION
%             %% UPDATE THE VALUES IN THE MAP
%             
%             
%         end        
%     end
    
    % -------------------------------------------------------------------------
    % Return key values
    % -------------------------------------------------------------------------
    function values = GetValues(obj, keys)
        % -------------------------------------------------------------------------
        % Prepare output
        % -------------------------------------------------------------------------         
        values = zeros(1, length(keys));
        
        % -------------------------------------------------------------------------
        % Find keys
        % -------------------------------------------------------------------------         
        validKeys = obj.data.isKey(num2cell(keys));
        if any(validKeys)
            values(validKeys) = cell2mat(obj.data.values(num2cell(keys(validKeys))));
        end
    end
    
    % -------------------------------------------------------------------------
    % Return Table
    % -------------------------------------------------------------------------
    function data = GetTable(obj)
        keys = cell2mat(obj.data.keys());
        values = cell2mat(obj.data.values());
        data = [keys; values;];
    end
    
    % -------------------------------------------------------------------------
    % Return keys
    % -------------------------------------------------------------------------
    function data = keys(obj)
        data = cell2mat(obj.data.keys());
    end

    % -------------------------------------------------------------------------
    % Return Values
    % -------------------------------------------------------------------------
    function data = values(obj)
        data = cell2mat(obj.data.values());
    end
    
    % -------------------------------------------------------------------------
    % Return length
    % -------------------------------------------------------------------------
    function numEntries = length(obj)
        numEntries = obj.data.Count;
    end
    
end
end