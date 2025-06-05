%% Demo cellpose integration


cellmask1 = imread(filePathToCellposeTif); % read in a tif file from cellpose
% here filePathToCellposeTif is the full filepath to the tif files created
% from cellpose. 

% if you have our current version of matlab-functions (from ORCA-public),
% you can Run's cellpose directly from Matlab and returns the results as a 
% matlab matrix using the following command:
%
% cellmask1 = RunCellpose(image2D); % 
% % -- requires you specify the filepaths to python and to cellpose 
% % this matlab function simply calls cellpose from your operating system
% % command line. 

% as a demo I'm going to use matlab's built in image of coins, but you
% should load your cellpose output files instead.  
sampleimage = imread('coins.png');
bwimage = imfill( sampleimage>100,'holes');
cellmask1 = bwlabel(bwimage);
figure(1); clf; imagesc(cellmask1 ); GetColorMap('distColorsW')

rna_xy = 255*rand(100,2); % the xy coordinates of your RNAs
% make sure these are in the same units as the image. ChrTracer for example
% converts DNA positions to nanometers, and you will need to convert back
% to pixels using the inverse of the conversion factor ChrTracer used in
% the first place. 

figure(1); hold on; plot(rna_xy(:,1),rna_xy(:,2),'k.');

[h,w] = size(cellmask1);
xy = round(rna_xy); % round to nearest pixel
% remove edge spots 
xy(xy(:,1)<1,:) =[];
xy(xy(:,2)<1,:) =[];
xy(xy(:,2)>h,:) =[];
xy(xy(:,1)>w,:) =[]; 

% covert x,y index to linear index
rna_idx = sub2ind([h,w], xy(:,2),xy(:,1));  % note that we need to put "row" (y-coord), before "column" (x-coodinate), so it's y,x not x,y)  

% get the cell ID for all RNA (DNA) spots
rna_cell_ID = cellmask1(rna_idx); 
% note, spots that fall outside the cell mask have ID=0.  

% an array that now has x-position, y-position, and rna_cell_ID
rna_table = [xy,rna_cell_ID] 


