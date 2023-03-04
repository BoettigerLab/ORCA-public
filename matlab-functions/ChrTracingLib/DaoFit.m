function data = DaoFit(frame,varargin)

% to use default parameters without editing, don't pass parameter arguments
% this should help with speed as well.

global daoSTORMexe defaultXmlFile scratchPath %#ok<GVMIS> 

defaultXmlFile = 'C:\Users\Alistair\Desktop\code\Archived\justSTORMfinder\STORMfinder\DefaultParameters\DaoFitPars_FirstFrame.xml';
defaultXmlFile = 'M:\2022-03-04_LiveTest48xTetO\b2_antifade_continuous\fitPars_test.xml'

% Parse default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'parsSaveFolder',{'scratchFolder','dataFolder',''},'scratchFolder'};
defaults(end+1,:) = {'parsFile','string',''};
defaults(end+1,:) = {'showPlots','boolean',true}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'threshold','nonnegative',14}; 
defaults(end+1,:) = {'iterations','nonnegative',4}; 
defaults(end+1,:) = {'background_sigma','nonnegative',8}; % 0 uses default in xml file  
defaults(end+1,:) = {'pixel_size','nonnegative',100}; % 0 uses default in xml file  
defaults(end+1,:) = {'start_frame','integer',-1}; % 0 uses default in xml file  
defaults(end+1,:) = {'max_frame','integer',1}; % 0 uses default in xml file  
defaults(end+1,:) = {'overwriteFit','boolean',false}; % 0 uses default in xml file  
% pars = ParseVariableArguments([],defaults);
pars = ParseVariableArguments(varargin,defaults,mfilename);

% to get all parameters
xmlPars = ReadXML(defaultXmlFile);
pars = JoinStructures(xmlPars.settings,pars,'conflict','keepSecond');

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
        xmlSavePath = [fileparts(frame),filesep];
    end
elseif length(frame)>2
    daxName = [scratchPath,'tempDax.dax'];
    WriteDax(frame,'saveFullName', daxName,'confirmOverwrite',false);
    xmlSavePath = scratchPath;
    pars.overwriteFit = true;
end

% ====== Create parameter file
% if no parameters passed
if ~isempty(pars.parsFile)
    fitPars = pars.parsFile;
else
    xmlPars = ReadXML(defaultXmlFile);
    updateFields = fieldnames(pars);
    new_values = cell(length(updateFields),1);
    for f=1:length(updateFields)
        xmlPars.settings.(updateFields{f}) = pars.(updateFields{f});
        new_values{f} =  pars.(updateFields{f});
    end
    scratchXmlFile = [xmlSavePath,'tempDaoPars.xml'];
    system(['copy ',defaultXmlFile, ' ',scratchXmlFile])
    UpdateTxtPars(scratchXmlFile,scratchXmlFile,updateFields,new_values,'verbose',false);  % 
    fitPars = scratchXmlFile; 
end


binName = regexprep(daxName,'.dax','.hdf5');
if exist(binName,'file')~=0
    if pars.overwriteFit
        delete(binName);
    else
        binName = IncrementSaveName(binName);   
    end
end

% call fit 
tic
cmdOut = [daoSTORMexe, ' --movie ',daxName,' --bin ', binName, ' --xml ', fitPars];
if pars.verbose
    disp(cmdOut);
    disp('now fitting image, may be slow for large images...');
end
system(cmdOut);
if pars.verbose
    tt =toc;
    disp(['fitting complete in ',num2str(tt),' seconds']);
end
%
nFrames = 1;
tic;
clear data;
data(nFrames).x = [];
data(nFrames).y = [];
data(nFrames).height = [];
data(nFrames).background = [];
data(nFrames).width = [];
for f=1:nFrames
    data(f).x = h5read(binName,['/fr_',num2str(f-1),'/x']);
    data(f).y = h5read(binName,['/fr_',num2str(f-1),'/y']);
    data(f).height = h5read(binName,['/fr_',num2str(f-1),'/height']);
    data(f).background = h5read(binName,['/fr_',num2str(f-1),'/background']);
    data(f).width = h5read(binName,['/fr_',num2str(f-1),'/xsigma']);
end
if pars.verbose
    tt =toc;
    disp(['data loaded in ',num2str(tt),' seconds']);
end


if pars.showPlots
    f=1;
    if haveDax
        daxFrame = ReadDax(frame,'startFrame',f,'endFrame',f,'verbose',false);
    else
        daxFrame = frame;
    end
    figure(1); clf; 
    imagesc(daxFrame); % IncreaseContrast(daxFrame,'high',.999)); 
    colorbar; colormap('gray');
    hold on; plot(data(f).x+1,data(f).y+1,'ro','MarkerSize',16);
end