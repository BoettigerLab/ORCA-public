function data = LoadFlyWig(datafile)
%% load HS-PolII data, see how it changes

% dataFolder = 'C:\Users\Alistair\Documents\Research\Projects\Chromatin\Data\2015-07-06-CorcesChIPseq\';
% wigFile = 'GSM1536013_RNAPII_HS_Rep1.wig';
% datafile = [dataFolder,wigFile];

% dm3 
fid = fopen(datafile);
fmt =['%f %f %*[^\n]'];
header = textscan(fid,'%s',1,'delimiter','\t');
polII_2L = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',2); 
polII_2Lhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_2R = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_2Rhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_3L = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_3Lhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_3R = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_3Rhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_4 = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_U = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_Uextra = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_X = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_Xhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_Yhet = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
polII_mito = textscan(fid,fmt,'CollectOutput',true,'delimiter','\t','HeaderLines',1); 
fclose(fid);
% figure(1); clf; plot(polII_2L{1}(:,1),polII_2L{1}(:,2),'b'); hold on;
% plot(polII_3L{1}(:,1),polII_3L{1}(:,2),'r');
data.chr2L.coords = polII_2L{1}(:,1);
data.chr2L.values = polII_2L{1}(:,2);
data.chr2R.coords = polII_2R{1}(:,1);
data.chr2R.values = polII_2R{1}(:,2);
data.chr3L.coords = polII_3L{1}(:,1);
data.chr3L.values = polII_3L{1}(:,2);
data.chr3R.coords = polII_3R{1}(:,1);
data.chr3R.values = polII_3R{1}(:,2);
data.chr4.coords = polII_4{1}(:,1);
data.chr4.values = polII_4{1}(:,2);
data.chrX.coords = polII_X{1}(:,1);
data.chrX.values = polII_X{1}(:,2);

