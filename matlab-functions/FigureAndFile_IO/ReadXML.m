function xmlStruct = ReadXML(xmlFile,varargin)
% Reads an xml file into a structure
% Nested xml data is converted into nested structure data
% 
% 
% xmlStruct = ReadXML('T:\ConvZscan_00.xml');
% xmlStruct = ReadXML('C:\Program Files\MATLAB\R2020b\toolbox\matlab\general\info.xml'); 
% 
% Example 1: repeat value of nested entry
% <parent>
% <child type="boolean">true</child>
% <child type="boolean">true</child>
% </parent>
% 
% data = ReadXML()
% struct data.parent with fields:
%       child: 1 (logical)
%     child_2: 1 (logical)

% xmlFile = 'T:\ConvZscan_00.xml';

pars.verbose = false;

% lazy way to read in a text file:
try
T = char(textread(xmlFile,'%s','delimiter','\n','whitespace','')); %#ok<DTXTRD>
catch er
    disp(er.message); 
    error(['could not read: ',xmlFile]);
end
[maxLines,~] = size(T);

% Parse the text file 
xmlStruct = struct();
try
    structLevel = {};
    nestCount = 1;
    f=1; % field counter
    for ln=1:maxLines
       [flagName,flagProps,flagStop] = ParseFlag(T,ln,pars);
       if ~isempty(flagName)
           flagName = regexprep(flagName,{'-'},{'_'}); % coerce valid field names
           flagOpen = ['<',flagName];
           flagClose = ['</',flagName '>'];
           isHeader = strcmp(flagOpen,'<?xml');
           isNestStop =  strcmp(flagName(1),'/');
           if ~isHeader && ~isNestStop
               [xClose,yClose] = FindFlagEnd(T,ln,flagClose,maxLines,pars);
               if ln~=yClose && xClose~=0
                   % start and end flags on different lines indicates branch structure
                   %------- check if flag exists, if so change name
                  newFlag = flagName;
                  try
                      cm = strcat( '.("', structLevel', '")' );
                      fieldExists = eval([ 'isfield(','xmlStruct',cat(2,cm{:}),',"',newFlag,'")']);
                  catch
                      fieldExists = false;
                  end
                  f = 1;
                  while fieldExists 
                      f=f+1;
                      newFlag = [flagName,'_',num2str(f)];
                      fieldExists = eval([ 'isfield(','xmlStruct',cat(2,cm{:}),',"',newFlag,'")']);
                  end
                  %-----------------
                   if ~isempty(flagProps)
                       if strcmp(flagProps{1,1},'validate') % record Validate status if specified
                           flagData = eval(lower( flagProps{1,2} )); %#ok<*NASGU>
                           cm = strcat( '.("', structLevel', '")' );
                           eval(['xmlStruct',cat(2,cm{:}),'.("',newFlag,'")','.("validate")=flagData;']); % a bit of a hack
                       end
                   else
                       cm = strcat( '.("', structLevel', '")' );
                       eval(['xmlStruct',cat(2,cm{:}),'.("',newFlag,'")','=struct();']); % a bit of a hack
                   end
                   % add to structLevel
                   structLevel{1,nestCount} = newFlag; %#ok<AGROW>
                   nestCount = nestCount + 1;
                   % now cycle on to the next line
               elseif ln==yClose && xClose~=0
                   % start and end flags on the same line indicates leaf level data
                  flagData =  T(ln,flagStop+1:xClose-1);
                  if ~isempty(flagProps)
                      if strcmp(flagProps{1,1},'type') % convert dataType if specified
                          if strcmp(flagProps{1,2},'boolean')
                            flagData = eval(lower( flagData )); %#ok<*NASGU>
                          elseif strcmp(flagProps{1,2},'int')
                              flagData = str2num(flagData); %#ok<ST2NM>
                          elseif strcmp(flagProps{1,2},'float')
                              flagData = str2num(flagData); %#ok<ST2NM>
                          end
                      end
                  end
                  %------- check if flag exists, if so change name
                  newFlag = flagName;
                  cm = strcat( '.("', structLevel', '")' );
                  try
                    fieldExists = eval([ 'isfield(','xmlStruct',cat(2,cm{:}),',"',newFlag,'")']);
                  catch
                      fieldExists = false;
                  end
                  f = 1;
                  while fieldExists 
                      f=f+1;
                      newFlag = [flagName,'_',num2str(f)];
                      fieldExists = eval([ 'isfield(','xmlStruct',cat(2,cm{:}),',"',newFlag,'")']);
                  end
                  %-----------------
                  cm = strcat( '.("', structLevel', '")' );
                  eval(['xmlStruct',cat(2,cm{:}),'.("',newFlag,'")','=flagData;']); % a bit of a hack
                  % disp(xmlStruct)
               elseif xClose==0
                   if pars.verbose
                       warning(['skipping ',flagName]);
                   end
               end
           end
           if isNestStop
               nestCount = nestCount-1;
               if ~isempty(structLevel)
                   structLevel(end) = [];
               end
           end
       end  
    end
catch er
    disp(er.message)
    disp('place debug here');
end

function [flagName,flagProps,flagStop] = ParseFlag(T,ln,pars) %#ok<INUSD>
    % parse xml flag
   flagStart = strfind(T(ln,:),'<');
   if length(flagStart) > 1 % take first occurrence if multiple
       flagStart = flagStart(1);
   end
   isComment = strcmp(T(ln,flagStart+1),'!');
   if isComment
       flagName = ''; flagProps = []; flagStop=0;
   else
       flagStop = strfind(T(ln,flagStart:end),'>') + flagStart-1;
       if length(flagStop) > 1 % take first occurrence if multiple
            flagStop = flagStop(1);
       end
       flagInfo = T(ln,flagStart+1:flagStop-1);
       flagParts = strsplit(flagInfo);
       % parse parts of flag Info
       flagName = flagParts{1}; % the first item is always the name
       numProps = length(flagParts)-1;
       flagProps = cell(numProps,2);
       for n=1:numProps
           typeValue = strsplit(flagParts{n+1},'=');
           flagProps{n,1} = typeValue{1};
           flagProps{n,2} = regexprep(typeValue{2},'"','');
       end
   end
   
function [xClose,yClose] = FindFlagEnd(T,ln,flagClose,maxLines,pars)
   yClose = ln;
   xClose = strfind(  T(yClose,:),flagClose);
   while isempty(xClose) && (yClose <= maxLines)     
       xClose = strfind(  T(yClose,:),flagClose);
       yClose = yClose+1;
   end
   if isempty(xClose)
       if pars.verbose
        warning(['Did not find end flag for ',flagClose]);
       end
       xClose = 0;
   end