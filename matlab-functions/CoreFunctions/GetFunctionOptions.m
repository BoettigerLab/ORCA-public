function [parStruct,parParts] = GetFunctionOptions(mFileName)
% Reads a matlab function formatted for ParseVariableArguments
%  returns the default parameter values
% 
% example application: create a parameter GUI for mFileName
% 
% parStruct = GetFunctionOptions('MyFxn.m');
% newPars = SimpleParameterGUI(parStruct);
% MyFxn('parameters',newPars);

% mFileName = 'ReadDax.m';
T = char(textread(mFileName,'%s','delimiter','\n','whitespace','')); %#ok<DTXTRD>

nLines = size(T,1);
for l=1:nLines
    parStart = contains(T(l,:),'defaults = cell(0,3)');
    if parStart
        parameterStartLine = l;
        break;
    end
end


for l=1:nLines
    parStart = contains(T(l,:),'ParseVariableArguments(');
    if parStart
        parameterEndLine = l;
        break
    end
end

nParLines = parameterEndLine - parameterStartLine + 1;
parParts = struct();
parParts(nParLines).name = '';
parParts(nParLines).type = '';
parParts(nParLines).value = '';
parParts(nParLines).comment = '';
for l=1:nParLines
   s = parameterStartLine + l; 
   d1 = strfind(T(s,:),'{');
   d2 = strfind(T(s,:),'}');
   c = strfind(T(s,:),'%');
   if isempty(c)
       c = length(T(s,:));
   end
   % some functions take cell array lists of values as parameters
   %   defaults(end+1,:) = {'method',{'linear','spline'},'linear'}
   t1 = [];
   t2 = [];
   if length(d1)>1
       t1 = d1(2); 
       d1 = d1(1);
   end
   if length(d2)>1
       t2 = d2(1);
       d2 = d2(end);
   end
   if ~isempty(d1)
        if ~isempty(t1)
           parParts(l).name = T(s,d1+2:t1-3);
           parParts(l).type = T(s,t1+1:t2-1);
           parParts(l).value =T(s,t2+3:d2-2);
           parParts(l).comment =T(s,d2+3:c);
        else
            % if the type is not comma separated, the first two comma
            % separated values are the name and type
            defaultInput = T(s,d1+1:d2-1);
            cs = strfind(defaultInput,',');
           parParts(l).name = defaultInput(1:cs(1)-1); 
           parParts(l).type = defaultInput(cs(1)+1:cs(2)-1);
           parParts(l).value = defaultInput(cs(2)+1:end);
           parParts(l).comment = T(s,c+1:end);
        end
        parParts(l).name = regexprep(parParts(l).name,"'",""); % names should be simple strings
        % parParts(l).value = regexprep(parParts(l).value,"'",""); % names should be simple strings
   end
end
% remove empty parameters;
parNames = {parParts.name};
parParts(cellfun(@isempty,parNames)) = [];
parNames = {parParts.name};
parValues = {parParts.value};

% convert to data structure, convert non-string data types/ 
nPars = length(parNames);
parStruct = struct();
for n=1:nPars
    try
        parStruct.(parNames{n}) = eval(parValues{n});
    catch
        parStruct.(parNames{n}) = parValues{n};
    end
%  % not needed, strings have double characters. 
%     if any(contains({'string','char'},parParts(n).type )) || contains(parParts(n).type,',')
%         parStruct.(parNames{n}) = parValues{n};
%     else
%         try
%             parStruct.(parNames{n}) = eval(parValues{n});
%         catch
%             parStruct.(parNames{n}) = parValues{n};
%         end
%     end
    
end