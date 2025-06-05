function [pntArray1b,pntArray2b] = ProcessTraces(daxFile1,varargin)
% 
% Output
% pntArray Matrix nTraces x nFrames x 6 data-values:
%    data-types = x_pix, y_pix, z_nm, t_frame, brightness, background 
%    Outputs two point arrays     
% 

% Defaults
npp = 108; % nm per pixel

defaults = cell(0,3);
% data loading
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'alignment_file','string',''}; % 
defaults(end+1,:) = {'chrom_correct_file','string',''}; % 
defaults(end+1,:) = {'bin_tag','string','_2d_iters'}; % 
defaults(end+1,:) = {'nFrames','integer',0}; % 0 for load all frames
% Trace identification 
defaults(end+1,:) = {'seedThresh','float',0.996}; % autoselect threshold for grouping localizations in z. 
defaults(end+1,:) = {'seedBinResolution','positive',2}; % autoselect threshold downsampling for grouping localizations in z. 
defaults(end+1,:) = {'doubletRadius','float',6}; % seed points within this distance from one another will be excluded to avoid hoping 
% matching
defaults(end+1,:) = {'maxSep','positive',25}; % max separation for spots to be considered a pair (in pixels)
defaults(end+1,:) = {'maxStep','positive',6}; % max step since last observed frame to consider linking (in pixels)
defaults(end+1,:) = {'maxDistToSeed','positive',16}; % max distance from seed point to link in a trace (in pixels)
% smoothing
defaults(end+1,:) = {'smooth','boolean',true}; % use smoothing?
defaults(end+1,:) = {'maxGap','integer',10}; %  max gap to attempt linear fill on (in frames)
defaults(end+1,:) = {'smoothWindow','integer',4}; %  moving average smoothing window size (in frames)
% z-fitting
defaults(end+1,:) = {'zFrames','integer',3}; % 
defaults(end+1,:) = {'zscan_nm','array',[-1000,0,1000]}; % 
defaults(end+1,:) = {'smoothZ','boolean',false}; % use smoothing?
% saveData 
defaults(end+1,:) = {'saveFolder','string',''}; %  save folder
defaults(end+1,:) = {'saveTables','boolean',true}; %  save folder
defaults(end+1,:) = {'saveMovies','boolean',true}; %  save folder
defaults(end+1,:) = {'cropRadius','integer',15}; % Radius of image around located spot to crop for movies.
pars = ParseVariableArguments(varargin,defaults,mfilename);



%% Main function
try

% === Get source folder and calibraion files
if isempty(pars.alignment_file)
    folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';
    pars.alignment_file = [folderCal2,'alignmentData.txt'];
end
if isempty(pars.chrom_correct_file)
    folderCal2 = 'M:\2023-08-26_130kb_rad21_dTag\';
    pars.chrom_correct_file = [folderCal2,'tform3D.mat'];
end
alignT = readtable(pars.alignment_file);
alignS = table2struct(alignT);

[dataFolder,daxName] = fileparts(daxFile1);
% create a save folder if none was passed
if isempty(pars.saveFolder)
    pars.saveFolder = [dataFolder,filesep,'Analysis\'];
    if ~exist(pars.saveFolder,'dir')
        mkdir(pars.saveFolder)
    end
end

daxName1 = [dataFolder,filesep,daxName,'.dax'];
daxName2 = regexprep(daxName1,'C1','C2'); %  [dataFolder,'36mW_0001_C2.dax'];
binName1 = regexprep(daxName1,'.dax',[pars.bin_tag,'.hdf5']); % [dataFolder,'36mW_0001_C1_2d_iters.hdf5'];  % 2d no iters
binName2 = regexprep(daxName2,'.dax',[pars.bin_tag,'.hdf5']);

% === load the images 
f = 1; 
[im1,info1] = ReadDax(daxName1,'startFrame',f,'endFrame',f+3,'verbose',false);  % note the +3 since its a 4 step series 
[im2,info2] = ReadDax(daxName2,'startFrame',f,'endFrame',f+3,'verbose',false);
im1 = max(im1,[],3);
im2 = max(im2,[],3);
im2 = ApplyReg(fliplr(im2),alignS);
im_f = cat(3,im1,im2); % higher contrast 
im_f = IncreaseContrast(im_f,'high',.9999,'low',.1);
figure(1); clf; Ncolor(im_f);  axis image;
[h,w] = size(im1);

% === Load the fit data
if pars.nFrames <= 0 
    nFrames = info1.number_of_frames;  %
else
    nFrames = pars.nFrames;
end
% nFrames = 200; % shorten for testing
fits1 = LoadHD5Fits(binName1,'nFrames',nFrames);    
fits2 = LoadHD5Fits(binName2,'nFrames',nFrames);

% === Apply camera alignment and chromatic correction
fits2c = Register2CamFits(fits2,'alignment_file',pars.alignment_file,'chrom_correct_file',pars.chrom_correct_file);

% === auto-select centers for dancing spots 
xy1 = TraceCenters({fits1.x},{fits1.y},'autoSelectThreshold',pars.seedThresh,'binResolution',pars.seedBinResolution);
xy2 = TraceCenters({fits2c.x},{fits2c.y},'autoSelectThreshold',pars.seedThresh,'binResolution',pars.seedBinResolution);

% === eliminate doublets
% challenge: averaged over all frames, the doublet won't have two centers. 
%  can't use just the first frame or first several, as not all spots are 
%  visible in frame 1:100, but may still be doublet)
dMap = squareform(pdist((xy1)));
dMap(dMap==0) = nan; % ignore self
nearestDot = min(dMap);  % 
xy1D = xy1(nearestDot<pars.doubletRadius,:);
xy1(nearestDot<pars.doubletRadius,:) = [];% eliminate both pairs
dMap = squareform(pdist((xy2)));
dMap(dMap==0) = nan; % ignore self
nearestDot = min(dMap);  % 
xy2D = xy2(nearestDot<pars.doubletRadius,:);
xy2(nearestDot<pars.doubletRadius,:) = [];% eliminate both pairs

% === pair selected spots
% keep only dancing spots that have partners in the other chromatic channel
[matched,cost] = MatchPoints(xy1,xy2,'maxDist',pars.maxSep);
m =matched(cost<pars.maxSep,:);      
nMatched = size(m,1);
xy1p = xy1(m(:,1),1:2);
xy2p = xy2(m(:,2),1:2); % indexed in order of starting point
    
% ---- view
figure(3); clf; Ncolor(im_f);  axis image;
hold on; plot(xy1p(:,1),xy1p(:,2),'y+');
hold on; plot(xy2p(:,1),xy2p(:,2),'gs');
sNum = cellstr(num2str((1:nMatched)'));
text(xy1p(:,1),xy1p(:,2),sNum,'color','r'); hold on;
% --------

% === Link spots into traces 
pntArray1 = LinkLiveStruct(fits1(1:nFrames),'seedPoints',xy1p,'maxStep',pars.maxStep,'maxDistToSeed',pars.maxDistToSeed); 
pntArray2 = LinkLiveStruct(fits2c(1:nFrames),'seedPoints',xy2p,'maxStep',pars.maxStep,'maxDistToSeed',pars.maxDistToSeed); 



f4 = figure(4); clf;
Ncolor(im_f);  axis image; hold on;
plot(cat(1,fits1.x),cat(1,fits1.y),'y+'); hold on;
plot(cat(1,fits2c.x),cat(1,fits2c.y),'gs');
hold on; plot(pntArray1(:,:,1)',pntArray1(:,:,2)','m.-');
hold on; plot(pntArray2(:,:,1)',pntArray2(:,:,2)','c.-');
hold on; plot(xy1D(:,1),xy1D(:,2),'mx');
hold on; plot(xy2D(:,1),xy2D(:,2),'cx');
hold on; plot(pntArray2(:,:,1)',pntArray2(:,:,2)','c.-');
text(xy1p(:,1),xy1p(:,2),sNum,'color','r'); hold on;
text(xy2p(:,1),xy2p(:,2),sNum,'color','c'); hold on;
% xlim([20,80]); ylim([1450,1500]);
% xlim([1340,1480]); ylim([260,390]);
% xlim([1,2024]); ylim([1,2024])
savefig(f4,[pars.saveFolder,'FOVoverview.fig'],'compact');
xlim([1340,1480]); ylim([260,390]); % zoom in for fun

figure(10); clf;  s=1:10;
subplot(2,1,1); plot(  (pntArray1(s,:,1) - nanmean(pntArray1(s,:,1),2))'  , '.-'); ylim([-16,16])
subplot(2,1,2); plot((pntArray2(s,:,1) - nanmean(pntArray2(s,:,1),2))', '.-' ); hold on; ylim([-16,16])
 
% a quick look at the pnt chain data 
f11 = figure(11); clf; ca = [-10,10];
subplot(1,4,1); imagesc(  squeeze(pntArray1(:,:,1)) -nanmean(squeeze(pntArray1(:,:,1)),2) );  colorbar; clim(ca); title('chain1 x')
subplot(1,4,2); imagesc(  squeeze(pntArray1(:,:,2)) -nanmean(squeeze(pntArray1(:,:,2)),2) ); colorbar; clim(ca);  title('chain1 y')
subplot(1,4,3); imagesc(  squeeze(pntArray2(:,:,1)) -nanmean(squeeze(pntArray2(:,:,1)),2) );  colorbar; clim(ca); title('chain2 x')
subplot(1,4,4); imagesc(  squeeze(pntArray2(:,:,2)) -nanmean(squeeze(pntArray2(:,:,2)),2) ); colorbar; clim(ca);  title('chain2 y')
GetColorMap('RedWhiteBlueSat');
savefig(f11,[pars.saveFolder,'TraceOverview.fig'],'compact');
% === Pair / match traces 
%  previous step matched based on all spots. This step re-matches the data
%  using only the linked spots. It is also necessary because the Linking
%  alogrithm is not guarenteed to preserve order [though if it doesn't drop
%  unlinked starting points, I think I could force it to prserve order].
p1c = squeeze(nanmean(pntArray1(:,:,1:2),2)); % average over frames
p2c = squeeze(nanmean(pntArray2(:,:,1:2),2)); % average over frames
p1c(isnan(p1c)) = 0;
p2c(isnan(p2c)) = 0;
[matched,cost] = MatchPoints(p1c,p2c,'maxDist',pars.maxSep);
m = matched(cost<pars.maxSep,:);      
nMatched = size(m,1);
pntArray1a = pntArray1(m(:,1),:,:);
pntArray2a = pntArray2(m(:,2),:,:);

subplot(1,2,2); Ncolor(im_f);  axis image;
hold on; plot(pntArray1a(:,:,1),pntArray1a(:,:,2),'y+');
hold on; plot(pntArray2a(:,:,1),pntArray2a(:,:,2),'gs');

% === noise filter
% Reduce time resolution to improve inference.  
% maxGap = 10;
% smoothWindow = 4; % 3 would integrates up and down 1 point, so the time resolution is effectively cut in 1/2. 
                 % 4 matches the z-scan duration of this dataset 
pntArray1b = pntArray1a;
pntArray2b = pntArray2a;
if pars.smooth
    for p=1:size(pntArray1a,1) 
        for d =1:2
            x1 = squeeze(pntArray1a(p,:,d));
            xx1 = fillmissing(x1,'linear','maxGap',pars.maxGap);
            xs1 = smooth(x1,pars.smoothWindow,'moving');
            xs1(isnan(xx1)) = nan;
            pntArray1b(p,:,d) = xs1;
    
            x1 = squeeze(pntArray2a(p,:,d));
            xx1 = fillmissing(x1,'linear','maxGap',pars.maxGap);
            xs1 = smooth(x1,pars.smoothWindow,'moving');
            xs1(isnan(xx1)) = nan;
            pntArray2b(p,:,d) = xs1;
        end
    end
end

% === re-order spots by data depth
isGood = ~isnan(pntArray1b(:,:,1) - pntArray2b(:,:,1));
[nGood,idx] = sort(sum(isGood,2),'descend');
pntArray1b = pntArray1b(idx,:,:);
pntArray2b = pntArray2b(idx,:,:);
% remove empty traces
pntArray1b(nGood<100,:,:) = [];
pntArray2b(nGood<100,:,:) = [];


% ====== Linked and matched
f4 = figure(4); clf;
Ncolor(im_f);  axis image; hold on;
plot(cat(1,fits1.x),cat(1,fits1.y),'y+'); hold on;
plot(cat(1,fits2c.x),cat(1,fits2c.y),'gs');
hold on; plot(pntArray1b(:,:,1)',pntArray1b(:,:,2)','m.-');
hold on; plot(pntArray2b(:,:,1)',pntArray2b(:,:,2)','c.-');
text(xy1p(:,1),xy1p(:,2),sNum,'color','r'); hold on;
text(xy2p(:,1),xy2p(:,2),sNum,'color','c'); hold on;
% xlim([20,80]); ylim([1450,1500]);
% xlim([1340,1480]); ylim([260,390]);
% xlim([1,2024]); ylim([1,2024])
savefig(f4,[pars.saveFolder,'FOV_Link_Overview.fig'],'compact');
xlim([1340,1480]); ylim([260,390]); % zoom in for fun

% ==== z calc
offFile =  regexprep(daxName1,'_C1.dax','.off');
offTable = ReadTableFile(offFile,'delimiter',' ');
zMat_1b = FitZTrace(pntArray1b,offTable,'parameters',pars);
pntArray1b(:,:,3) = zMat_1b;
zMat_2b = FitZTrace(pntArray2b,offTable,'parameters',pars);
pntArray2b(:,:,3) = zMat_2b;



% ==== save table
disp(['writing data tables'])
if pars.saveTables
    [nCells,nFrames,nMeasurements] = size(pntArray1b);
    for c=1:nCells
        x1_nm =  pntArray1b(c,:,1)'*npp;
        y1_nm =  pntArray1b(c,:,2)'*npp;
        z1_nm =  pntArray1b(c,:,3)'; % already in nm
        h1 =  pntArray1b(c,:,5)';
        bkd1 = pntArray1b(c,:,6)';
        t1 = (1:length(x1))';
        x2_nm =  pntArray2b(c,:,1)'*npp;
        y2_nm =  pntArray2b(c,:,2)'*npp;
        z2_nm =  pntArray2b(c,:,3)';
        h2 =  pntArray2b(c,:,5)';
        bkd2 = pntArray2b(c,:,6)';
        t2 = pntArray2b(c,:,6)';
        distXY_nm = sqrt(sum(  (pntArray1b(c,:,1:2) - pntArray2b(c,:,1:2)).^2, 3))'.*npp;
        ctable = table(x1_nm,y1_nm,z1_nm,h1,bkd1,t1,x2_nm,y2_nm,z2_nm,h2,bkd2,t2,distXY_nm);
        writetable(ctable,[pars.saveFolder,'traj_',num2str(c,'%30d'),'.txt']);
    end
end

% === save images


%% assemble movies
% load the images

if pars.saveMovies
    disp(['loading raw movie data...'])
    nS = nCells;
    xf = nan(nS,1);
    yf = nan(nS,1);
    r = pars.cropRadius;
    nT = ceil(nFrames/pars.zFrames);
    im5D = zeros(2*r+1,2*r+1,pars.zFrames,2,nT,nS,'uint16');
    t=0;
    figure(1); clf; 
    k = 1; % counter of each frame in a batch
    b = 1; % counter for the batch number
    batchSize = 100;
   for f = 1:1:nFrames
            % read movie in chunks.
            %   when running in parpool this should reduce the number of 
            %   simultaneous read/write commands trying to pull from the 
            %   same disk.  
            if k==1
                ff = batchSize*(b-1)+k; 
                fe = min(ff+batchSize-1,nFrames);
                % disp(ff)
                im1_dax = ReadDax(daxName1,'startFrame',ff,'endFrame',fe,'verbose',false); 
                im2_dax = ReadDax(daxName2,'startFrame',ff,'endFrame',fe,'verbose',false);
                if pars.verbose
                    disp([num2str(ff/nFrames),'% movie loaded'])
                end
            else
                if k==batchSize
                    k=0; % reset
                    b=b+1;
                end
            end
             % disp([f,k,b,ff])
            im1 = im1_dax(:,:,k+1);
            im2 = im2_dax(:,:,k+1);
%             im1 = ReadDax(daxName1,'startFrame',f,'endFrame',f,'verbose',false); 
%             im2 = ReadDax(daxName2,'startFrame',f,'endFrame',f,'verbose',false);
            im2 = ApplyReg(fliplr(im2),alignS);
            im3 = IncreaseContrast(cat(3,im1,im2),'high',.99995,'low',.001); % constant contrast
            [h,w] = size(im1);
        o = rem(f,pars.zFrames); 
        if o==0
            o=pars.zFrames; % keeping indexes straight
        elseif o==1
            t=t+1; % starting new series
        end
        for s=1:nS % s =4
            try
           % update centering
                x = round(pntArray1b(s,f,1)); % in pixels :)
                y = round(pntArray1b(s,f,2));
                if ~isnan(x)
                    xf(s) = x;
                    yf(s) = y;
                end
            % plot
            if ~isnan(xf(s))
                x1 = max([1,xf(s)-r]);
                x2 = min([w,xf(s)+r]);
                y1 = max([1,yf(s)-r]);
                y2 = min([h,yf(s)+r]);
                imS = im3(y1:y2,x1:x2,:);
                xL = x2-x1+1;
                yL = y2-y1+1;
                if (x2 > 2304) || (y2 > 2304)  % shouldn't be necessary, 
                    disp(([y1,y2,x1,x2]))
                    disp(([h,w]))
                    disp('debug here')
                end
                im5D(1:yL,1:xL,o,:,t,s) = imS; 
            end
            catch er
                warning(['error processing spot ',num2str(s), ' frame ',num2str(f)])
                disp(([y1,y2,x1,x2]))
                warning(er.message);
                warning(er.getReport);
                disp('debug here')
            end
        end
        % set(gcf,'color','k');
        % pause(.01);
        k=k+1;
    end
    
    
    %% save files
    saveMovieFolder = [pars.saveFolder,'movies\'];
    disp(['saving cropped movies as i5d...'])
    if ~exist(saveMovieFolder,'dir')
        mkdir(saveMovieFolder);
    end
    
    for s=1:nS
        im5Dout = im5D(:,:,:,:,:,s);
        % [dH,dW,dZ,dC,dT] = size(im5Dout);
         WriteImage5D(im5Dout,[saveMovieFolder,'spot_',num2str(s,'%03d'),'.i5d'],'lociX',xf(s),'lociY',yf(s));
    end

end


%% example on plotting movie

%  %% Load and View Saved Files
%  s=4;
% [im5,im5_info] = ReadImage5D([saveFolder,'spot_',num2str(s,'%03d'),'.i5d']);
% for t=1:dT  % t=1
%     figure(1); clf; 
%     for z=1:pars.zFrames
%         subplot(1,pars.zFrames,z);
%         im = squeeze(im5(:,:,z,:,t));
%         Ncolor(im); title(t); 
%         hold on;
%         xf(s) = im5_info.lociX;
%         yf(s) = im5_info.lociY;
%         % some troubleshooting plots
%       plot([fits1(f).x]-xf(s)+r+2,[fits1(f).y]-yf(s)+r+2,'ys'); hold on; % all spots 1
%     plot([fits2c(f).x]-xf(s)+r+2,[fits2c(f).y]-yf(s)+r+2,'bs'); hold on; % all spots 2  (to verify linking)
%     plot(pntArray1a(:,f,1)-xf(s)+r+2,pntArray1a(:,f,2)-yf(s)+r+2,'yo'); hold on; % linked spots from all arrays (for cross talk checking)
%     % the key plots 
%     plot(pntArray1a(s,f,1)-xf(s)+r+2,pntArray1a(s,f,2)-yf(s)+r+2,'w+'); hold on; % linked spots 1
%     plot(pntArray2a(s,f,1)-xf(s)+r+2,pntArray2a(s,f,2)-yf(s)+r+2,'b+'); hold on; % linked spots 2
%         title(['f',num2str(f),'  s',num2str(s),'  rb',num2str(pntArray1a(s,f,5),4),'  gb',num2str(pntArray2a(s,f,5),4)],'color','w')
%     end
%     pause(.01); 
% end
catch er
    warning(er.message)
    warning(er.getReport)
    disp('place debug here')
end

