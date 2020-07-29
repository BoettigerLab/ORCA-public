I = imread('ellipse_wiki.png');
E = edge(rgb2gray(I),'canny');

%% override some default parameters
params.minMajorAxis = 20;
params.maxMajorAxis = 50;

% note that the edge (or gradient) image is used
bestFits = ellipseDetection(mask1, params);

fprintf('Output %d best fits.\n', size(bestFits,1));

figure(4); clf;
image(mask1); hold on;
%ellipse drawing implementation: http://www.mathworks.com/matlabcentral/fileexchange/289 
a=1;
ellipse(bestFits(a,3),bestFits(a,4),bestFits(a,5)*pi/180,bestFits(a,1),bestFits(a,2),'r');
