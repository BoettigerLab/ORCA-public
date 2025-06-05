function output = AlignColorsStepGUI(gui_step,pars,data,varargin)
% 
% 
% defaults = cell(0,3);
% parameters = ParseVariableArguments(varargin,defaults,mfilename);

%%  Basic Step GUI executable
%  The idea of this code organization is to easily wrap a step-by-step
%  block organized script into a simple step-by-step parameter adjusting
%  GUI.  
%  The function takes a step number, a task ('task') 
% 
%  
% 
%% Align Colors Step GUI
% data is taken with both lasers on with multi-color (tetraspec) beads
%  the images in C2 should be primarily the blue colors (bounced by the long pass dichroic) 
%  the images in C1 should be primarily the red colors (passed by the dichroic)
%  consequently these need both a chormatic alignment and a xy-theta
%  registration. "camera registration"
% 
% The images are taken in 3D. 
%  the current default is:
%  bead movies are taken using software z-scan, with +/- 4 um scan range
%  and 5 dead frames at the start (which should get dropped). 
%  - these properties should be read from the XML and Offset file directly
%  
%
%% example
% general    stpGUI1 = MiniStepGUI();
% this version  clrApp1 = AlignColorsStepApp();
% 
%% some functions of interest
%  DaoFit.m
% LoadHD5Fits.m
% LoadDax
%
%

output = []; 

%% Step 1: set data path, save path, load the data and display,
if gui_step == 1
    % folderCal2 = 'I:\Jude\2023-10-10_400kb_test\';
    if isempty(pars) % if no parameters are passed get them for this step
        pars.beadC1_dax = '';
        pars.directions = 'Set data path, save path, load the data and display';
        output = pars;
   
    else
        
        if isempty(pars.beadC1_dax)
            [file,folder] = uigetfile('*.dax','select a bead C1 dax image');
            beadC1_dax = [folder,file];
            pars.dataPath = folder; %        'I:\Jude\2023-10-12_400kb_test\' ;  % INPUT
        else
            beadC1_dax = pars.beadC1_dax;
            pars.dataPath = [fileparts(beadC1_dax),filesep]; % (fileparts does not add the filesep, but uigetfile does  
        end
        beadC2_dax = regexprep(beadC1_dax,'C1','C2');
        beads1 = ReadDax(beadC1_dax);
        beads2 = ReadDax(beadC2_dax);

        % read xml file with scope parameters 
        xmlFile = regexprep(beadC1_dax,'.dax','.xml');
        xmlFile = regexprep(xmlFile,'_C1','');
        xmlStruct = ReadXML(xmlFile);

        % read offset file 
        offsetFile = regexprep(beadC1_dax,'.dax','.off');
        offsetFile = regexprep(offsetFile,'_C1','');
        offsetTable = ReadTableFile(offsetFile,'delimiter',' ');

        % display max project 
        d1 = max(beads1,[],3);
        d2 = max(beads2,[],3);
        imD = cat(3,d1,fliplr(d2));
        figure(2); clf; Ncolor(IncreaseContrast(imD,'high',.9995,'low',.2)); 

        output.dataFolder = pars.dataPath; 
        output.beadC1_dax = beadC1_dax;
        output.beadC2_dax = beadC2_dax;
        output.dax1 = beads1;
        output.dax2 = beads2;  % these are small enough it's worth storing just 2 complete z-scans 
        output.dax_xml = xmlStruct; 
        output.dax_offsetTable = offsetTable;
        output.rawOverlay = imD;
    end
end


%% Align the cameras
if gui_step == 2
    if isempty(pars)
        pars.directions = 'Compute camera registration';
        pars.alignmentFileName = 'alignmentData.txt'; % 
        pars.angles = 0.50:.05:0.80; % 
        pars.fineCenter = [1000,1000]; % 
        output = pars;
    else
        alignment_file = [data.dataFolder,pars.alignmentFileName];
        if exist(alignment_file,'file')
            alignT = readtable(alignment_file);
            alignC = table2struct(alignT);
        else   
            f30 = figure(30); clf;
            alignC3 = CorrAlignFast(data.rawOverlay(:,:,1),data.rawOverlay(:,:,2),...
                'angles',pars.angles,'maxSize',300,'fineCenter',pars.fineCenter);  %  fast align first
            [h,w,~] = size(data.rawOverlay);
            alignC = alignC3;
            alignC.h = h;
            alignC.w = w;
            alignC.flip = {'flipLR_none_none'}; % testing putting all these in one 
            alignT = struct2table(alignC);
            %  Save the data!
            writetable(alignT,alignment_file);
            save([data.dataFolder,'alignment.mat'],"alignC");
            SetFigureSavePath(data.dataFolder); 
            SaveFigure(f30,'name','alignmentIm','formats',{'fig','png'});
        end
        
        % note, fliplr has already been applied to rawOverlay
        cc = ApplyReg(data.rawOverlay(:,:,2),alignC,'invert',false);
        imC = cat(3,data.rawOverlay(:,:,1),cc);
        figure(2); clf; Ncolor(IncreaseContrast(imC,'high',.999,'low',.2)); colorbar;
        xlim([900,1100]); ylim([900,1100]);

        % save ouputs
        output.alignC = alignC; 
        output.dataOverlay = imC;
    end
end



%% test how well we are doing with registration of other images 

if gui_step == 3
    if isempty(pars)
        pars.directions = 'Validate registration on another image.  Note, this step is optional,select Next Step to skip.';
         [file,folder] = uigetfile([data.dataFolder,'*.dax'],'select a control bead C1 dax image');
        if ~isempty(file)
            pars.dataPath = folder; %        'I:\Jude\2023-10-12_400kb_test\' ;  % INPUT
            pars.beadDaxFiled = file; % 'beads_0001_C1.dax';
        end
        output = pars;
    else
        beadC1 = [pars.dataPath,pars.beadDaxFiled];
        beadC2 = regexprep(beadC1,'C1','C2');
        data_C1 = ReadDax(beadC1);
        data_C2 = ReadDax(beadC2);
        d1 = max(data_C1,[],3);
        d2 = max(data_C2,[],3);
        dd = ApplyReg(fliplr(d2),data.alignC,'invert',false);
        imD = cat(3,d1,dd);
        figure(2); clf; Ncolor(IncreaseContrast(imD,'high',.9995,'low',.2)); 
        output.controlOverlay = imD; 
    end
end


%% beads are beautiful, lets fit!
if gui_step == 4
% Step 1, fit the dax files and load the data
% this is slow. 
    if isempty(pars)
        pars.directions = 'compute fits. warning, this step is slow, wait till windows close';
        pars.beadC1_dax = data.beadC1_dax;
        pars.beadC2_dax = data.beadC2_dax;
        pars.binC1 = '';
        pars.binC2 = '';
        pars.parsFile = ''; % 'M:\2023-08-26_130kb_rad21_dTag\fitPars_allframe_iters_C2.xml';
        pars.overwrite = 'resume';
        output = pars;
    else
        if isempty(pars.parsFile)
            pars.parsFile = [data.dataFolder,'fitPars_beads.xml'];
            if ~exist(pars.parsFile,'file')
                [file,folder] = uigetfile('*.xml','select an xml DaoSTORM parameters file');
                pars.parsFile =[folder,file];
            end
        end
        if isempty(pars.binC1)
            output.binC1 = regexprep(pars.beadC1_dax,'.dax','.hdf5');
        else
            output.binC1 = pars.binC1;
        end
        if isempty(pars.binC2)
            output.binC2 = regexprep(pars.beadC2_dax,'.dax','.hdf5');
        else
            output.binC2 = pars.binC2;
        end
        fits_c1 = [];
        fits_c2 = [];
        if ~exist(output.binC1,'file') || strcmp(pars.overwrite,'overwrite')
            fits_c1 = DaoFit(pars.beadC1_dax,'parsFile',pars.parsFile,'runExternal',true,'overwrite',pars.overwrite);
        else
            disp('existing data found, proceed to next step')
        end
        if ~exist(output.binC2,'file') || strcmp(pars.overwrite,'overwrite')
            fits_c2 = DaoFit(pars.beadC2_dax,'parsFile',pars.parsFile,'runExternal',true,'overwrite',pars.overwrite);
        end
        output.fits_c1 = fits_c1;
        output.fits_c2 = fits_c2; 
    end
end

%% validate 2D fits
if gui_step == 5
    if isempty(pars)
        pars.directions = 'load and plot fits';
        output = pars;
    else
         fits_c1 = LoadHD5Fits(data.binC1);
         fits_c2 = LoadHD5Fits(data.binC2);

         % note, the c2 channels has an extra-dead frame at the start.
         fits_c2([fits_c2.frame]==1) = [];
         ExpandAll = @(temp) temp{:}; % need a helper function
         [fits_c2.frame] = ExpandAll(num2cell([fits_c2.frame]-1));
        
    
        % (this requires the camera registration computed above)
        fits_c2a = Register2CamFits(fits_c2,'alignment_file',data.alignC);
        
        % figure(2); clf; 
        % subplot(1,2,1); plot(cat(1,fits_c1.x)',cat(1,fits_c1.y)','b+'); hold on;
        % subplot(1,2,2); plot(cat(1,fits_c1.x)',cat(1,fits_c1.height)','b+'); hold on;
        % subplot(1,2,1); plot(cat(1,fits_c2a.x)',cat(1,fits_c2a.y)','ro'); 
        % xlabel('x (pixels)'); ylabel('y (pixels)');
        % subplot(1,2,2); plot(cat(1,fits_c2a.x)',cat(1,fits_c2a.height)','ro');
        % xlabel('x (pixels)'); ylabel('brightness of fit (gaussian height)')
        
        
        dd = ApplyReg(data.rawOverlay(:,:,2),data.alignC,'invert',false);
        imD = cat(3,data.rawOverlay(:,:,1),dd);
        figure(1); clf; Ncolor(IncreaseContrast(imD,'high',.9995,'low',.2)); hold on;
        plot(cat(1,fits_c1.x)',cat(1,fits_c1.y)','b+'); hold on;
        plot(cat(1,fits_c2a.x)',cat(1,fits_c2a.y)','ro'); 
        xlabel('x (pixels)'); ylabel('y (pixels)');
        title('fit overlay (zoom in to validate spots)'); 

        % select variables to save
        output.fits_c1 = fits_c1;
        output.fits_c2 = fits_c2;
        output.fits_c2a = fits_c2a;
        data.dataOverlay = imD;
    end
end
%% FitZ
if gui_step == 6
    if isempty(pars)
        pars.directions= 'Fit z-data';
        pars.useOffsets = true;
        pars.removePoorFits = true;
        pars.maxZ  = 10e3; % z values outside this range will be removed
        pars.minZ  = 0; % z values outside this range will be removed
        pars.overwrite = false; 
        try
            pars.startZ = data.dax_xml.settings.focuslock.software_z_scan.deadtime + 1;
            pars.nppZ = data.dax_xml.settings.focuslock.software_z_scan.step_size;
            pars.zrange = data.dax_xml.settings.focuslock.software_z_scan.range;
            pars.endZ = pars.startZ + 2*pars.zrange/pars.nppZ;
        catch
            pars.startZ = 6; % 5 dead frames
            pars.endZ = 86; % 80 data frames
        end
        output = pars; 
    else
        % inputs (translate) 
        fits_c1 = data.fits_c1;
        fits_c2a = data.fits_c2a;
        % remove start and end beads, in case it is confusing to the
        % matching algorithm to have the same beads twice  (?)
        fits_c1([fits_c1.frame]<pars.startZ) = [];
        fits_c1([fits_c1.frame]>pars.endZ) = [];
        fits_c2a([fits_c2a.frame]<pars.startZ) = [];
        fits_c2a([fits_c2a.frame]>pars.endZ) = [];

        % search for existing fit data
        ztable1_name = regexprep(data.beadC1_dax,'.dax','_zfits.csv');
        ztable2_name = regexprep(data.beadC2_dax,'.dax','_zfits.csv');
        if exist(ztable1_name,'file' ) && exist(ztable2_name,'file') && ~pars.overwrite 
            disp('found existing data, loading these');
            fits_c1z = readtable(ztable1_name);
            fits_c2z = readtable(ztable2_name);
        else

            
    
            % deal with truncated offset table (Hal error, why does this happen?)  
            o = data.dax_offsetTable.offset;
            if pars.endZ > length(o)
                pars.endZ = length(o);
            end
            if pars.useOffsets
                o = data.dax_offsetTable.offset;
                o(o==0) = nan;
                zz = pars.startZ:pars.endZ; % in units of frames 
                xx = o(zz); % in units of offset 
                zz(isnan(xx)) = nan;
                figure(1); clf;  % figure(2); clf;
                yyaxis left; plot(xx,zz,'.'); xlabel('offset'); ylabel('frame');
                goodData = ~isnan(xx);

                 pp = polyfit(xx(goodData),zz(goodData),1);
                 zf = pp(2)*xx+pp(1); % best fit (in units of offset)
                 hold on; plot(xx,zf,'-');  ylim([0,length(o)]);
                 z_per_frame = zf*pars.nppZ;
                 z_per_frame = fillmissing(z_per_frame,'linear');
                 yyaxis right; plot(xx,z_per_frame,'>');    ylim([0,length(o)*pars.nppZ]);
                ylabel('z (nm)')
                title('offset to z conversion');
                z_per_frame_tot = [zeros(pars.startZ-1,1); z_per_frame; zeros(length(o)-pars.endZ,1)];
            else
                z_per_frame_tot = (0:pars.nppZ:(length(o)-1)*pars.nppZ)'; % ignore offsets 
            end      
            data.dax_offsetTable.z_nm = z_per_frame_tot; 
            writetable(data.dax_offsetTable,[data.dataFolder,'offsets_nm.txt']);
            % 
         
            disp('computing z-fits...'); % these can also be slow, maybe we should save files?  ideally we'd write this into the .hdf5 file, but I haven't hacked into taht one yet.  
            tic
            fits_c1z = DaoFitZ(fits_c1,'sigma',pars.nppZ*5,'seedThresh',.05,'z_per_frame',z_per_frame,'output','table'); % ,'multiFit',true
            fits_c2z = DaoFitZ(fits_c2a,'sigma',pars.nppZ*5,'seedThresh',.05,'z_per_frame',z_per_frame,'output','table'); % ,'multiFit',true
            toc
            % remove nans
            noData = isnan(fits_c1z.z);
            fits_c1z(noData,:) = []; 
            noData = isnan(fits_c2z.z);
            fits_c2z(noData,:) = []; 
            
               % === save the z-fits
            writetable(fits_c1z,ztable1_name); 
            writetable(fits_c2z,ztable2_name);      
        end

     % remove poor fits
        if pars.removePoorFits
            figure(1); clf; subplot(2,1,1); hist(fits_c1z.z,0:100:20e3);
            poorFit = fits_c1z.ciratio >1  | fits_c1z.z > pars.maxZ  | fits_c1z.z < pars.minZ ;
            fits_c1z(poorFit,:) = []; % 
           figure(1); subplot(2,1,2);  hist(fits_c1z.z,0:100:20e3);
            poorFit = fits_c2z.ciratio>1  | fits_c2z.z > pars.maxZ | fits_c2z.z < pars.minZ;
            fits_c2z(poorFit,:) = []; % toss the bad fits;
        end

        % (mostly useful for troubleshooting)
        figure(1); clf; 
        plot(cat(1,fits_c1.x)',cat(1,fits_c1.y)','k+','MarkerSize',3); hold on;
        plot(cat(1,fits_c2a.x)',cat(1,fits_c2a.y)','ro','MarkerSize',3); 
         plot(fits_c1z.x,fits_c1z.y,'g>'); hold on; 
         plot(fits_c2z.x,fits_c2z.y,'m<'); 
         xlabel('x (pixels)'); xlabel('y (pixels)');
         title('compare Z-fit seedpoints with data');
         legend('chn 1 spots','chn 2 spots','chn 1 centers','chn 2 centers');

        % % for troubleshooting 
        % data = stpGUI1.data;     pars = stpGUI1.pars.step6; 
        % fits_c1z = data.fits_c1z; fits_c2z = data.fits_c2z;
        
        figure(3); clf; % let's show this as an overlay instead
        subplot(1,2,1); plot(fits_c1z.x,fits_c1z.y,'b+'); hold on; 
        subplot(1,2,2); plot(fits_c1z.x,fits_c1z.z,'b+'); hold on;
        subplot(1,2,1); plot(fits_c2z.x,fits_c2z.y,'ro'); 
        xlabel('x (pixels)'); ylabel('y (pixels)');
        subplot(1,2,2); plot([fits_c2z.x]',[fits_c2z.z]','ro');
        xlabel('x (pixels)'); ylabel('z (nm)');
        % 
        % figure(4); clf; % this is best for troubleshooting
        % plot(log10([fits_c1z.h])',[fits_c1z.z]','b+'); hold on;
        % plot(log10([fits_c2z.h])',[fits_c2z.z]','ro'); hold on;
        % xlabel('fit height'); ylabel('z position')    
        
    % display 4D overlay for validation
        try
            xs = 500:700; % 1000:1200;
            ys = 500:700;
            zs = pars.startZ:pars.endZ; % 100;
  
            dax2lr = fliplr(data.dax2);
            data.dax2a = dax2lr;
            for z=1:size(dax2lr,3)
                data.dax2a(:,:,z) = ApplyReg(dax2lr(:,:,z),data.alignC,'invert',false);
            end
    
            x1 = fits_c1z.x;
            y1  = fits_c1z.y;
            x2 = fits_c2z.x;
            y2  = fits_c2z.y;
            % fz1 = ([fits_c1z.z])/pars.nppZ - pars.startZ+1;  %  these don't quite line up
            % fz2 = ([fits_c2z.z])/pars.nppZ - pars.startZ+1;  %  these don't quite line up
            fz1 =(fits_c1z.z)/pars.nppZ;  %  these don't quite line up
            fz2 = (fits_c2z.z)/pars.nppZ;  %  these don't quite line up
            inBox1 = x1 > xs(1) & x1<xs(end) & y1 > ys(1) & y1 < ys(end); %  & fz1 > zs(1) & fz1 < zs(end);
            inBox2 = x2 > xs(1) & x2<xs(end) & y2 > ys(1) & y2 < ys(end); %  & fz2 > zs(1) & fz2 < zs(end);
            
            im4D = cat(4,data.dax1(ys,xs,zs),data.dax2a(ys,xs,zs));
            figure(1); clf; ProjectIm4D(im4D);
            subplot(1,3,1);  hold on;  plot(x1(inBox1)-xs(1)+1,y1(inBox1)-ys(1)+1,'ro','MarkerSize',18);
            subplot(1,3,1);  hold on;  plot(x2(inBox2)-xs(1)+1,y2(inBox2)-ys(1)+1,'co','MarkerSize',18);
            subplot(1,3,2);  hold on;  plot(y1(inBox1)-ys(1)+1,fz1(inBox1)-zs(1)+1,'ro','MarkerSize',18);
            subplot(1,3,2);  hold on;  plot(y2(inBox2)-ys(1)+1,fz2(inBox2)-zs(1)+1,'co','MarkerSize',18);
            subplot(1,3,3);  hold on;  plot(x1(inBox1)-xs(1)+1,fz1(inBox1)-zs(1)+1,'ro','MarkerSize',18);
            subplot(1,3,3);  hold on;  plot(x2(inBox2)-xs(1)+1,fz2(inBox2)-zs(1)+1,'co','MarkerSize',18);

        catch er
            disp(er.message);
            disp('can not display 4D overlay');
        end

        dd = ApplyReg(data.rawOverlay(:,:,2),data.alignC,'invert',false);
        imD = cat(3,data.rawOverlay(:,:,1),dd);
        figure(2); clf; Ncolor(IncreaseContrast(imD,'high',.9995,'low',.2)); hold on;
        plot(cat(1,fits_c1.x)',cat(1,fits_c1.y)','b+'); hold on;
        plot(cat(1,fits_c2a.x)',cat(1,fits_c2a.y)','ro'); 

        % outputs (translate)
        output.fits_c1z = fits_c1z;
        output.fits_c2z = fits_c2z;
        disp('step completed')
    end
end
%% Matchpoints 3d
% note, z is already in nm from the DaoFitZ, but xy are still in pixels
if gui_step == 7
    if isempty(pars)
        pars.directions= 'Match points in 3D';
        pars.maxDist = 800; % in nm
        try
            obj1Data = strsplit(data.dax_xml.settings.mosaic.obj1,',');
            nppXY = 1e3*str2double(obj1Data{2});
            pars.nppXY = nppXY; %  108;
        catch
            pars.nppXY = 108;
        end
        output = pars; 
    else
        disp('matching points...')
        % inputs (translate)
        fits_c1z = data.fits_c1z;
        fits_c2z = data.fits_c2z;

        % parameters (translate)
        maxDist = pars.maxDist; % input parameter

        tic
        % p1 = [ [fits_c1z.x]'*pars.nppXY,[fits_c1z.y]'*pars.nppXY,[fits_c1z.z]'];
        % p2 = [ [fits_c2z.x]'*pars.nppXY,[fits_c2z.y]'*pars.nppXY,[fits_c2z.z]'];
        p1 = [ fits_c1z.x*pars.nppXY,fits_c1z.y*pars.nppXY,fits_c1z.z];
        p2 = [ fits_c2z.x*pars.nppXY,fits_c2z.y*pars.nppXY,fits_c2z.z];
        [matched,cost] = MatchPoints(p1,p2,'maxDist',maxDist);
        m =matched(cost<maxDist,:); 
        pm1 = p1(m(:,1),:);
        pm2 = p2(m(:,2),:);
        toc
        
        figure(3); clf; 
        subplot(1,2,1); plot(pm1(:,1),pm1(:,2),'k+'); hold on; 
        subplot(1,2,2); plot(pm1(:,1),pm1(:,3),'k+'); hold on;
        subplot(1,2,1); plot(pm2(:,1),pm2(:,2),'ro'); 
        xlabel('x (pixels)'); ylabel('y (pixels)');
        subplot(1,2,2); plot(pm2(:,1),pm2(:,3),'ro');
        xlabel('x (pixels)'); ylabel('z (nm)');
        
        figure(4); clf; 
        plot3(p1(:,1),p1(:,2),p1(:,3),'b.'); hold on;
        plot3(pm1(:,1),pm1(:,2),pm1(:,3),'k+'); hold on;
        plot3(pm2(:,1),pm2(:,2),pm2(:,3),'ro'); hold on;
        plot3([pm1(:,1),pm2(:,1)]',[pm1(:,2),pm2(:,2)]',[pm1(:,3),pm2(:,3)]','-'); hold on;
        xlabel('x (pixels)'); ylabel('y (pixels)'); zlabel('z (nm)');

% xl = get(gca,'xlim');  yl = get(gca,'ylim');  xlim(xl); ylim(yl);
        % ouptuts (translate)
        output.pm1 = pm1; % points matched chn1. 
        output.pm2 = pm2; % points matched chn2
        disp('step completed')

    end
end
%% apply 
if gui_step == 8

    if isempty(pars)
        pars.directions= 'Compute chromatic correction maps.';
        pars.maxDist = 600; 
        output = pars; 
    else
        disp('computing chromatic maps...')
        % inputs (translate)
        % nppXY = pars.nppXY;
        % nppZ = pars.nppZ; 
        maxDist = pars.maxDist;
    
        % pm1 = data.pm1;
        % pm2 = data.pm2;
        % pts1 = [pm1(:,1:2)*nppXY,pm1(:,3)*nppZ];
        % pts2 = [pm2(:,1:2)*nppXY,pm2(:,3)*nppZ];
        tic
        tform3D = Polymap3D(data.pm1,data.pm2,'max2D',maxDist,'max3D',10*maxDist,'polyOrder',2,'bins',0:20:1000,'units','nm');
        toc
        
        % % save
%         tformPars.npp = pars.nppXY; % xy
        % tformPars.nppZ = pars.nppZ; % xy
        output.tform3D = tform3D;
%        output.tformPars = tformPars;
        disp('step completed');
    end
end

%% save
if gui_step == 9
    if isempty(pars)
        pars.directions = 'save results';
        pars.saveName = 'chromatic';
        output = pars;
    else
        tform3D = data.tform3D; 
    save([data.dataFolder,'tform3D','_',pars.saveName,'.mat'],"tform3D");
    % save(['F:\DefaultPars\','tform3D.mat'],"tformPars","tform3D");
    f1 = figure(23);
        SetFigureSavePath(data.dataFolder);
    SaveFigure(f1,'name',['tform_results','_',pars.saveName],'formats',{'png'});
    f2 = figure(22);
    SaveFigure(f2,'name',['tform_map','_',pars.saveName],'formats',{'png'})
    end
end
%% DONE!
% 
% some testing now
