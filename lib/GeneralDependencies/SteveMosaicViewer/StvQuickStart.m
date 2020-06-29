function [xu,yu,xp,yp,z,mag,mtiles,M] = StvQuickStart(mosaicPath,varargin)
% converts .stv files to .mat files and saves a StvQuickStart file of the
% position data.  


defaults = cell(0,3);
defaults(end+1,:) = {'restart','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% mosaicPath = 'E:\Alistair\2017-06-19_Calibration\561_Mosaic\561_mosaic.msc'
[mosaicFolder,mosaicName] = fileparts(mosaicPath);


% Read in .mat files and get list
%-----------------------------------
qstart = dir([mosaicFolder,'quickstart.mat']);

if isempty(qstart) || pars.restart
    mtiles = dir([mosaicFolder,filesep,mosaicName,'*.mat']);
    M = length(mtiles); 
    if M==0
        disp(['found ',num2str(M),' tiles in folder ',mosaicFolder]); 
        disp(['Attempting to convert .stv files to .mat']); 
        ConvertStv2Mat(mosaicPath);
        mtiles = dir([mosaicFolder,filesep,mosaicName,'*.mat']);
        M = length(mtiles); 
    end
    if M==0
        error('no .mat mosaic files found');
    end
    if pars.verbose
        disp(['found ',num2str(M),' tiles in folder ',mosaicFolder]);
    end

    xu = zeros(M,1);
    yu = zeros(M,1);
    xp = zeros(M,1);
    yp = zeros(M,1); 
    z = zeros(M,1); 
    mag = zeros(M,1);

    for m=1:M
         load([mosaicFolder,filesep,mtiles(m).name]);
             xp(m) = x_pix;
             xu(m) = x_um;
             yp(m) = y_pix;
             yu(m) = y_um;
             mag(m) = magnification;
             z(m) = zvalue;
             if pars.verbose
                disp(['sorting data... ',num2str(100*m/M,3),'%']);
             end
    end

    save([mosaicFolder,'quickstart.mat'],'xu','yu','xp','yp','z','mag','mtiles','M'); 
    if pars.verbose
        disp('writing quikstart.mat file');
    end
else
    load([mosaicFolder,'quickstart.mat']); 
end
    