function infoFile = WriteDax(dax,varargin)
% infoFile = WriteDax(dax)
% Writes the matrix "dax" to a .dax file
% automatically completes necessary information in the .inf file
%
%% Example
% WriteDax(dax,'folder',saveFolder,'daxName','myImage');
% 
%% Default Parameters
% defaults(end+1,:) = {'folder', 'string', scratchPath}; % 
% defaults(end+1,:) = {'daxName', 'string', 'temp'}; % (leave off .dax)
% defaults(end+1,:) = {'verbose', 'boolean', true}; % 
% 
global scratchPath

defaults = cell(0,3);
defaults(end+1,:) = {'folder', 'string', scratchPath}; % 
defaults(end+1,:) = {'daxName', 'string', 'temp'}; % (leave off .dax)
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'dataType','string','little endian'};
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% enforce no 'dax' in the inf name
parameters.daxName = regexprep(parameters.daxName,'.dax','');

[yDim,xDim,nFrames] = size(dax); 

                 infoFile.localName= [parameters.daxName,'.inf'];
                 infoFile.localPath= parameters.folder;
                  infoFile.uniqueID= 0; % 7.3616e+05;
                      infoFile.file= '';
              infoFile.machine_name= 'matlab-storm';
           infoFile.parameters_file= '';
             infoFile.shutters_file= '';
                  infoFile.CCD_mode= 'frame-transfer';
                 infoFile.data_type= ['16 bit integers (binary, ',parameters.dataType,')'];
          infoFile.frame_dimensions= [xDim,yDim]; %[yDim xDim];
                   infoFile.binning= [1 1];
                infoFile.frame_size= xDim*yDim;
    infoFile.horizontal_shift_speed= 0;
      infoFile.vertical_shift_speed= 0;
                infoFile.EMCCD_Gain= 0;
               infoFile.Preamp_Gain= 0;
             infoFile.Exposure_Time= 0;
         infoFile.Frames_Per_Second= 0;
        infoFile.camera_temperature= 0;
          infoFile.number_of_frames= nFrames;
               infoFile.camera_head= 'matlab';
                    infoFile.hstart= 1;
                      infoFile.hend= xDim;
                    infoFile.vstart= 1;
                      infoFile.vend= yDim;
                 infoFile.ADChannel= 0;
                   infoFile.Stage_X= 0;
                   infoFile.Stage_Y= 0;
                   infoFile.Stage_Z= 0;
               infoFile.Lock_Target= 0;
                  infoFile.scalemax= 1844;
                  infoFile.scalemin= 100;
                     infoFile.notes= '';
                     
 WriteDAXFiles(dax,infoFile,'verbose',parameters.verbose);
 