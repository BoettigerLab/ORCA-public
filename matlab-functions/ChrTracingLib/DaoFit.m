function data = DaoFit(frame,varargin)
% 
% to use default parameters without editing, don't pass parameter arguments
% this should help with speed as well.

global daoSTORMexe defaultXmlFile scratchPath matlabFunctionsPath %#ok<GVMIS> 
global condaPrompt
if isempty(daoSTORMexe)
    warning('global variable daoSTORMexe not defined')
    error('check DaoSTORM is installed and added to your matlab startup and filepath')
end
if isempty(scratchPath)
    warning('global variable scratchPath not defined')
    error('Please define a scratchPath to avoid creating stray daxfiles')
end

% clean start
tempFiles = {[scratchPath,'tempDax.dax'],...
    [scratchPath,'tempDax.hdf5'],...
    [scratchPath,'tempDax.csv'],...
    [scratchPath,'tempDax.inf'],...
     [scratchPath,'tempDaoPars.xml']};
for f=1:length(tempFiles)
    if exist(tempFiles{f},'file')
        delete(tempFiles{f});
    end
end

% defaultXmlFile = 'C:\Users\Alistair\Desktop\code\Archived\justSTORMfinder\STORMfinder\DefaultParameters\DaoFitPars_FirstFrame.xml';
% defaultXmlFile = 'M:\2022-03-04_LiveTest48xTetO\b2_antifade_continuous\fitPars_test.xml'
% defaultXmlFile = [matlabFunctionsPath,'ChrTracingLib\fitPars_test.xml'];

% Parse default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'parsSaveFolder',{'scratchFolder','dataFolder',''},'dataFolder'};
defaults(end+1,:) = {'binName','string',''}; % leave empty to autopopulate. otherwise enter a full filepath here  
defaults(end+1,:) = {'parsFile','string',''}; %  use this parameter file for analysis
defaults(end+1,:) = {'modParsFile','string',''};  % use this parameter file for analysis, but update 
defaults(end+1,:) = {'updateParsStruct','struct',[]}; % structure of parameters to update
defaults(end+1,:) = {'showPlots','nonnegative',1}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'threshold','nonnegative',14}; 
defaults(end+1,:) = {'iterations','nonnegative',20}; 
defaults(end+1,:) = {'sigma','nonnegative',1}; % PSF sigma 
defaults(end+1,:) = {'background_sigma','nonnegative',8.0}; % 0 uses default in xml file
defaults(end+1,:) = {'foreground_sigma','nonnegative',2.0}; % 0 uses default in xml file
defaults(end+1,:) = {'camera_offset','nonnegative',100}; % 0 uses default in xml file  
defaults(end+1,:) = {'find_max_radius','nonnegative',5}; % 0 uses default in xml file  
defaults(end+1,:) = {'model',{'2d','2dfixed','3d','Z'},'2d'}; % 0 uses default in xml file  
defaults(end+1,:) = {'pixel_size','nonnegative',100}; % 0 uses default in xml file  
defaults(end+1,:) = {'start_frame','integer',-1}; % 0 uses default in xml file  
defaults(end+1,:) = {'max_frame','integer',1}; %  
defaults(end+1,:) = {'x_start','integer',0}; %  
defaults(end+1,:) = {'x_stop','integer',2304}; %  
defaults(end+1,:) = {'y_start','integer',0}; %   
defaults(end+1,:) = {'y_stop','integer',2304}; % 
defaults(end+1,:) = {'overwrite',{'resume','increment','overwrite'},'resume'};
defaults(end+1,:) = {'method',{'h5read','csv'},'h5read'};   % h5read (needed for DaoFitZ, which is currently written only to use the crazy   
% defaults(end+1,:) = {'overwriteFit','boolean',false}; % 0 uses default in xml file 
% defaults(end+1,:) = {'rerun','boolean',false}; % 0 uses default in xml file 
defaults(end+1,:) = {'runExternal','boolean',false}; % 0 uses default in xml file 
defaults(end+1,:) = {'Hidden','boolean',false}; % 0 uses default in xml file 
% pars = ParseVariableArguments([],defaults);
pars = ParseVariableArguments(varargin,defaults,mfilename);

% to get all parameters
if isempty(pars.parsFile) && isempty(pars.modParsFile)
    currXml = defaultXmlFile; 
    xmlPars = ReadXML(defaultXmlFile);
    pars = JoinStructures(xmlPars.settings,pars,'conflict','keepSecond');
elseif ~isempty(pars.modParsFile)
    currXml = pars.modParsFile; 
    xmlPars = ReadXML(currXml);
    if isempty(pars.updateParsStruct)
        updatePars = struct();
        for j=1:2:length(varargin)
            updatePars.(varargin{j}) = varargin{j+1};
        end
    else
        updatePars = pars.updateParsStruct;
    end
    parsTemp = JoinStructures(xmlPars.settings,updatePars,'conflict','keepSecond'); % update xml with any passed options
    pars = JoinStructures(pars,parsTemp,'conflict','keepSecond'); % also keep the matlab specific pars in the structure
else
    currXml = pars.parsFile; 
    xmlPars = ReadXML(currXml);
    pars = JoinStructures(xmlPars.settings,pars,'conflict','keepFirst');
end

%%


% ===== Test if argument was a Dax File or an image
haveDax = false;
% frame = 'M:\2022-03-04_LiveTest48xTetO\b2_antifade_continuous\48xTetO_B2_antifade_0001.dax';
% frame =
% 'K:\2022-04-21_Repeat_TwistRNA_CTCF_Rad_E14\Hyb_107\maxProjections\Max750_00.dax'; % myc exon   
% frame = 'K:\2022-04-21_Repeat_TwistRNA_CTCF_Rad_E14\Hyb_109\maxProjections\Max750_00.dax'; % nanog exon 
% frame = 'K:\2022-04-21_Repeat_TwistRNA_CTCF_Rad_E14\Hyb_105\maxProjections\Max750_01.dax'; % Fgf5 exon 
% frame = 'K:\2022-04-21_Repeat_TwistRNA_CTCF_Rad_E14\Hyb_098\maxProjections\Max750_01.dax'; % Prdm14 exon 
if isstr(frame)
    if strcmp(frame(end-3:end),'.dax')
        haveDax = true;
        daxName = frame;
        if strcmp(pars.parsSaveFolder,'dataFolder')
            xmlSavePath = [fileparts(frame),filesep];
        elseif strcmp(pars.parsSaveFolder,'scratchPath')
            xmlSavePath = scratchPath;
        end
    end
elseif length(frame)>2
    daxName = [scratchPath,'tempDax.dax'];
    WriteDax(frame,'saveFullName', daxName,'confirmOverwrite',false,'verbose',false);
    xmlSavePath = scratchPath;
    pars.overwriteFit = true;
    pars.rerun = true;
end

% ====== Create parameter file
% if no parameters passed
if ~isempty(pars.parsFile)
    fitPars = pars.parsFile;
elseif ~isempty(pars.modParsFile)  && ~isempty(pars.updateParsStruct)
    xmlPars = ReadXML(pars.modParsFile);
    updateFields = fieldnames(pars.updateParsStruct);
    new_values = cell(length(updateFields),1);
    for f=1:length(updateFields)
        xmlPars.settings.(updateFields{f}) = pars.updateParsStruct.(updateFields{f});
        new_values{f} =  pars.updateParsStruct.(updateFields{f});
    end
    % ---- this part could be a helper function
    scratchXmlFile = [xmlSavePath,'tempDaoPars.xml'];
        system(['copy "',currXml, '" "',scratchXmlFile,'"']);
    UpdateTxtPars(defaultXmlFile,scratchXmlFile,updateFields,new_values,'verbose',false);  % 
    fitPars = scratchXmlFile; 
    % --------
else
    xmlPars = ReadXML(currXml);
    updateFields = fieldnames(pars);
    new_values = cell(length(updateFields),1);
    for f=1:length(updateFields)
        xmlPars.settings.(updateFields{f}) = pars.(updateFields{f});
        new_values{f} =  pars.(updateFields{f});
    end
    % -----could be helper function
    scratchXmlFile = [xmlSavePath,'tempDaoPars.xml'];
    system(['copy "',currXml, '" "',scratchXmlFile,'"']);
    UpdateTxtPars(defaultXmlFile,scratchXmlFile,updateFields,new_values,'verbose',false);  % 
    fitPars = scratchXmlFile; 
    % -------
end



% 3 options:
%     (1) detect existing file and resume if analysis not complete
%     (2) create new fit (hdf5) field using increment name
%     (3) overwrite 
%  
if isempty(pars.binName)
    binName = regexprep(daxName,'.dax','.hdf5');
    if exist(binName,'file')~=0
        if strcmp(pars.overwrite,'overwrite')
            delete(binName);
        elseif strcmp(pars.overwrite,'resume')
            disp('existing file found, will load this')
        else  % pars.overwrite 'increment'
            binName = IncrementSaveName(binName);   
        end
    end
else
    binName = pars.binName;
end

% call fit 
% cmdOut = [condaPrompt ,daoSTORMexe, ' --movie ',daxName,' --bin ', binName, ' --xml ', fitPars];
cmdOut = [daoSTORMexe, ' --movie ',daxName,' --bin ', binName, ' --xml ', fitPars];
if pars.verbose
    disp(cmdOut);
    disp('now fitting image, may be slow for large images...');
end
if pars.runExternal
    if pars.verbose
        disp('running externally');
    end
    data = SystemRun(cmdOut,'Hidden',pars.Hidden);  % actually returns proc. 

elseif ~pars.runExternal
    tic
    system(cmdOut);
    if pars.verbose
        tt =toc;
        disp(['fitting complete in ',num2str(tt),' seconds']);
    end
    % If not running external, we can display data when complete
    if nargout > 0 ||  pars.showPlots
        data = LoadHD5Fits(binName,'parameters',pars);
    end
    if pars.verbose
        tt =toc;
        disp(['data loaded in ',num2str(tt),' seconds']);
    end
    if pars.showPlots
        f=1;
        ff = pars.start_frame + f -1;
        if pars.max_frame > pars.start_frame
            fe = pars.max_frame;
        else
            fe = ff;
        end
        if haveDax
            daxFrame = ReadDax(frame,'startFrame',ff,'endFrame',fe,'verbose',false);
        else
            daxFrame = frame;
        end
        daxFrame = max(daxFrame,[],3);
        figure(pars.showPlots); clf; 
        imagesc(daxFrame); % IncreaseContrast(daxFrame,'high',.999)); 
        colorbar; colormap('gray');
        xx = cat(1,data.x) + 1;
        yy = cat(1,data.y) + 1;
        hold on; plot(xx,yy,'ro','MarkerSize',16);
    end
end