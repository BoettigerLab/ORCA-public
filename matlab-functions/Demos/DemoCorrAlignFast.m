
% Demo CorrAlignFast
close all;

folder = fullfile(matlabroot, '\toolbox\images\imdata');
im = imread([folder,filesep,'rice.png']);
im1 = IncreaseContrast(im,'low',.7,'high',.9998);

figure(1); clf; imagesc(im1);



av.xshift = 23.5;
av.yshift = -12.5;
av.theta = -2;
av.rescale = .9;
im2 = ScaleRotateShift(im1,av,'invert',true);
im3 = ScaleRotateShift(im1,av);


figure(3); clf; 
subplot(1,2,1); Ncolor(cat(3,im1,im2));
subplot(1,2,2); Ncolor(cat(3,im1,im3));

scales = .8:.05:1;
angles = -2:1;

tic
avOut = CorrAlignFast(im1,im2,... 
    'scales',scales,...
    'angles',angles,...
    'maxSize',100,...
    'showExtraPlot',true,...
    'fineUpsample',2,...
    'minFineImprovement',-inf)
toc   

tic
avOut2 = CorrAlignRotateScale(im1,im2,...  % very slow
    'scales',scales,...
    'angles',angles)
toc