function spotTable = CT4_BlankSpotTable(N)

% N = pars.nFits*nDataChns;
%=== Main table ================
spotID = nan(N,1);
traceID = nan(N,1);
x_um = nan(N,1);  % spot position relative to fid (useful for distance filtering)  
y_um = nan(N,1);  % spot position relative to fid
z_um = nan(N,1);  % spot position relative to fid
x_fid = nan(N,1); % fid position relative to FOV (useful for FOV specific errors like chromatic aberration and illumination)
y_fid = nan(N,1); % fid position relative to FOV
z_fid = nan(N,1); % fid position relative to FOV
barcodeID = nan(N,1);
cellID = nan(N,1);
%= additional spot data==========
fitQuality = nan(N,1); % false if fit outputs = inputs, denotes failed fit
x_shift_um = nan(N,1);
y_shift_um = nan(N,1);
z_shift_um = nan(N,1);
shift_score = nan(N,1);
spt_brightness = nan(N,1);
x_um_ci95 = nan(N,1);
y_um_ci95 = nan(N,1);
z_um_ci95 = nan(N,1);
spt_bkd = nan(N,1);
spt_bkd_xtilt = nan(N,1);
spt_bkd_ytilt = nan(N,1);
spt_sXY = nan(N,1);
spt_sZ = nan(N,1);
fid_brightness = nan(N,1);
dat_chn = nan(N,1);
% ---- data common to all spots
stageX = nan(N,1); % fov position relative to stage/sample
stageY = nan(N,1); % fov position relative to stage/sample
fov  = nan(N,1);
hyb = nan(N,1);

spotTable = table(spotID,traceID,x_um,y_um,z_um,barcodeID,cellID,...
        fitQuality,x_fid,y_fid,z_fid,x_shift_um,y_shift_um,z_shift_um,shift_score,spt_brightness,x_um_ci95,y_um_ci95,z_um_ci95,...
        spt_bkd,spt_bkd_xtilt,spt_bkd_ytilt,spt_sXY,spt_sZ,fid_brightness,dat_chn,...
        stageX,stageY,fov,hyb);