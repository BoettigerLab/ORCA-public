function pars = GetScopeSettings(infoStruct,varargin)
% pars = GetScopeSetting(infoFile);
% pars = GetScopeSetting(infoStruct);
% pars = GetScopeSetting([],'scope','scope1');
% pars = GetScopeSetting([],'parameters',parsIn);  where parsIn has field
% pars.scope1
% Note, the varargin input allows this function to update parameter values
% with just the scope name.  By keeping this in one place, if we change the
% Parameters for any given scope it is easy to update.  
% 
scopeList = {'scope1','scope2','scope3','scope2b','other','autoDetect'};

defaults = cell(0,3);
defaults(end+1,:) = {'scope',scopeList,'autoDetect'};  %
defaults(end+1,:) = {'transpose','boolean',false}; 
defaults(end+1,:) = {'fliplr','boolean',false}; 
defaults(end+1,:) = {'flipud','boolean',false}; 
defaults(end+1,:) = {'nmPixXY','positive',100}; 
defaults(end+1,:) = {'pix_to_mm','positive',10}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

usedAuto = false;
% try to guess scope format automatically from the InfoFile
if ~contains(pars.scope,scopeList)
    warning(['scope name ',pars.scope ' not recognized']);
    disp(' valid options are: ');
    disp(scopeList);
end
    

% scope1 info file
% infoFile = ReadInfoFile('N:\Derek\20200929_L10-b2-E14mESC_im6\DNA_Expt\Hyb_001\ConvZscan_00.dax')
% 
% % read from infoFile.settings.  This should be mor elegant for the future
% pars.scope = infoFile.settings.setup_name;
% pars.fliplr = infoFile.settings.mosaic.flip_horizontal;
% pars.flipud = infoFile.settings.mosaic.flip_vertical;
% pars.transpose = infoFile.settings.mosaic.transpose;
% pars.nmPixXY = 154;  % should hard record this in xml file for the future. The spotcounter has a nmPix field
% pars.pix_to_mm = 1/pars.nmPixXY*1e6; % = 152.67 ~= 154;  but pretty close!

if strcmp(pars.scope,'autoDetect')
    usedAuto = true;
    % try to read the InfoFile. Through an error if failed
    err = false;
    if ischar(infoStruct)
        try
            infoStruct = ReadInfoFile(infoStruct);
        catch
            err = true;
        end
    end
    if ~isstruct(infoStruct)
        err = true;
    end
    if err
        error('could not read infoStruct. Must be either the full file path to an .inf file associated with a dax movie or a matlab structure loaded from such a file');
    end

    % use size
    imHeight = infoStruct.frame_dimensions(2); % size(dax,1);
    imWidth = infoStruct.frame_dimensions(1); % size(dax,1);
    if imHeight == 1024
        pars.scope = 'scope1';
    elseif imHeight == 1536
        pars.scope = 'scope2';
    elseif imHeight == 2048
        pars.scope = 'scope3';
    elseif imHeight == 1540 
        pars.scope = 'scope2b';
    else
        pars.scope = 'other';
    end
end

if pars.verbose
    if usedAuto
        disp(['autoDetect determined the data is from system: ',pars.scope]);
    else
        disp(['System was instructed data is from: ',pars.scope]);
    end
end


switch pars.scope
    % Scope1= old values =T T F.  10/27/20 -> T F T - Mosaic Builder Valid 
    % Scope2= old values =T F T.  07/13/20 -> T T F
    % Scope3=T F T

    case 'scope1'              % T  - Mosaic Builder valid
        pars.transpose = true;% 
        pars.fliplr = false;   % 
        pars.flipud = true;   % 
        pars.nmPixXY = 154;
        pars.pix_to_mm = 6.55; % 1/pars.nmPixXY*1e6
        disp([pars.transpose pars.fliplr pars.flipud]);
    case 'scope2'
        pars.transpose = true;  % T  F  T    T     Settings=T (only good if Steve does it the same way)
        pars.fliplr = true;     % F  F  T    F         Settings=T
        pars.flipud = false;    % F  T  F    T          Settings=F
        pars.nmPixXY = 154;
        pars.pix_to_mm = 6.55;
    case 'scope2b'
        % Jeremy 12/11
        pars.transpose = true; 
        pars.fliplr = false;    
        pars.flipud = true;   
        pars.nmPixXY = 154;
        pars.pix_to_mm = 6.55;        
    case 'scope3'
        pars.transpose = true;
        pars.fliplr = false;
        pars.flipud = true;
        pars.nmPixXY = 108;
        pars.pix_to_mm = 9.3398; % 4.5935; %  6.55;
     case 'other'
        opars.transpose = true;
        opars.fliplr = true;
        opars.flipud = false;
        opars.nmPixXY = 154;
        opars.pix_to_mm = 6.55;
        newPars = SimpleParameterGUI(opars);
        pars.transpose = newPars.transpose;
        pars.fliplr = newPars.fliplr;
        pars.flipud = newPars.flipud;
        pars.nmPixXY = newPars.nmPixXY;
        pars.pix_to_mm = newPars.pix_to_mm;
    case 'otherwise'  
        warning('did not autodetect scope format. Using default or requested scope parameters.');
        warning('Erroneous display might result if these parameters are wrong.'); 
end
    
pars.imHeight = imHeight;
pars.imWidth = imWidth;
