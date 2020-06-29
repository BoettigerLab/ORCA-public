function [mosaicImage,mosaicPars,aveIm] = LoadStvTiles(varargin)
% [mosaicImage,mosaicPars] = LoadStvTiles() 
%  -- Loads the steve mosaic currently in global variable stvfile and plots
%  the results in figure(1) for the first 1000 frames.   
% [mosaicImage,mosaicPars] = LoadStvTiles(mosaicFile) 

%%
global stvfile 
% flattenPath = 'E:\Alistair\2017-06-19_Calibration\ave488.tif';
% flattenPath = 'E:\Alistair\2017-06-19_Calibration\ave561.tif';
% flattenPath = 'E:\Alistair\2017-06-19_Calibration\ave647.tif';

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'showall', 'boolean', true}; 
defaults(end+1,:) = {'showbox', 'boolean', false};  % plot box over seed position 
defaults(end+1,:) = {'showNumbers','boolean',true};
defaults(end+1,:) = {'shrk', 'positive', 1}; 
defaults(end+1,:) = {'position', 'array', [0,0]};  % seed position
defaults(end+1,:) = {'N', 'positive', 2000};  % max number of tiles around seed position 
defaults(end+1,:) = {'frameSize', 'positive', 1024}; 
defaults(end+1,:) = {'flattenPath','string',''};
defaults(end+1,:) = {'flatten','boolean',false};
defaults(end+1,:) = {'smoothFlatten','nonnegative',0};
defaults(end+1,:) = {'restart','boolean',false};
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------

if isempty(varargin)
    mosaicFile = stvfile;
    varInput = [];
else
    mosaicFile = varargin{1};
    varInput = [];
end
if length(varargin) > 1
    varInput = varargin(2:end);
end

parameters = ParseVariableArguments(varInput, defaults, mfilename);

% default parameters
position = parameters.position;

if ~isempty(parameters.flattenPath)
    aveIm = imread(parameters.flattenPath);
    % load(parameters.flattenPath,'aveIm');
    if parameters.verbose
       disp('loading background image for flat-field correction'); 
    end
else
    aveIm = [];
end


%% Main Function

[mosaicFolder,mosaicName] = fileparts(mosaicFile); %#ok<*ASGLU>
mosaicFolder = [regexprep(mosaicFolder,'\','/'),'/']; 

[xu,yu,xp,yp,z,mag,mtiles,M] = StvQuickStart(mosaicFile,'parameters',parameters);
%%
N = min(parameters.N,M); 
[frames,~] = knnsearch([xu,yu],position(1,:),'k',N);

% get dimensions of grid to plot:
if parameters.showall
    magmin = 1/min(mag(frames));
else
    magmin = 1;
end
    
xmin = min(xp(frames)-magmin*parameters.frameSize);
xmax = max(xp(frames)+magmin*parameters.frameSize);
ymin = min(yp(frames)-magmin*parameters.frameSize);
ymax = max(yp(frames)+magmin*parameters.frameSize);

xs = round(xmax - xmin); 
ys = round(ymax - ymin);
X = round(xp-xmin);
Y = round(yp-ymin);

shrk = parameters.shrk;
while max([ys,xs]) > 1E7
   warning('Mosaic is larger than 1E7 x 1E7.  Downsampling...');
    shrk = shrk*2; 
    xs = round(xs/shrk); 
    ys = round(ys/shrk); 
    X = round(X/shrk);
    Y = round(Y/shrk); 
    disp(['downsampled mosaic ' ,num2str(shrk), ' fold']); 
end

% This is heavy on memory for large N or for scaling up low res images.
% mosaicImage = zeros(ys,xs,1,'uint16');
mosaicImage = nan(ys,xs,1,'single');

allIm = cell(length(frames),1);
for k = 1:length(frames)
    i = frames(k);
    load([mosaicFolder,filesep,mtiles(i).name]);
    im = squeeze(data);
    allIm{k} = im;
    
    if parameters.verbose
        if rem(k,5) == 0
          disp(['loading data ', num2str(100*k/length(frames),3), '% done']);
        end
    end
end
allIm = cat(3,allIm{:});

if parameters.flatten
    % aveIm = nanmedian(allIm,3);
    aveIm = nanmedian(allIm,3);
    [h,w] = size(aveIm);
    if parameters.smoothFlatten > 0
        aveIm = imgaussfilt(aveIm,h/parameters.smoothFlatten);
    end
end



for k = 1:length(frames)
    i = frames(k);
    im = squeeze(allIm(:,:,k));
    % for plotting the 4x images, need to scale up:
    if parameters.showall
        if mag(i) ~= 1
            im = imresize(im,1/mag(i));
        end 
    end
    
 % rescale image if requested
     if shrk ~= 1
        im = imresize(im,1/shrk); 
     end
    [h,w,~] = size(im);
    h2 = floor(h/2); % avoid integer operands for colon operator warning
    w2 = floor(w/2);    
    y1 = Y(i)+1-h2;
    y2 = Y(i)+h2;
    x1 = X(i)+1-w2;
    x2 = X(i)+w2;
    
    if ~isempty(aveIm)
        im = 2^10*single(im)./ single(imresize(aveIm,[h,w]));
        % im = uint16(2^10* double(im)./double(imresize(aveIm,[h,w])));
        % disp('applying normalization');
    end
    
  try  
      mosaicImage(y1:y2,x1:x2) = nanmean( cat(3,mosaicImage(y1:y2,x1:x2), im'),3);
     % mosaicImage(y1:y2,x1:x2) = im';
  catch
      mosaicImage(y1:y2-1,x1:x2-1) = nanmean( cat(3,mosaicImage(y1:y2-1,x1:x2-1), im'),3);
      % mosaicImage(y1:y2-1,x1:x2-1) = im';
  end
  %   mosaicImage(x1:x2,y1:y2) = fliplr( flipud(im) );

    if parameters.verbose
        if rem(k,5) == 0
          disp(['rebuilding mosaic ', num2str(100*k/length(frames),3), '% done']);
        end
    end
end

%---- convert back to uint16
mosaicImage = uint16(mosaicImage);

%--------- show plot
clf;
imOut = imadjust(mosaicImage,stretchlim(mosaicImage,0.001));
imagesc(imOut); 
colormap(gray(2^8));
set(gcf,'color','w');
axis image;

%-------------------- compute pixel to um conversion    
% %figure(2); plot(xu(mag==1),xp(mag==1),'k.');
% mx = ( max(xp(mag==1)) - min(xp(mag==1)) )/( max(xu(mag==1)) - min(xu(mag==1)) ) ;
% my = ( max(yp(mag==1)) - min(yp(mag==1)) )/( max(yu(mag==1)) - min(yu(mag==1)) ) ;

mx = ( max(xp) - min(xp) )/( max(xu) - min(xu) ) ;
my = ( max(yp) - min(yp) )/( max(yu) - min(yu) ) ;

mosaicPars.mx = mx; 
mosaicPars.my = my;
mosaicPars.xmin = xmin;
mosaicPars.ymin = ymin;
mosaicPars.frameSize = parameters.frameSize/mag(i);

box_coords = PositionToBox(position,mosaicPars,'frameSize',parameters.frameSize*1/mag(i));
numPos = size(box_coords,1);

if parameters.showbox
    hold on;
    for n=1:numPos
      rectangle('Position',box_coords(n,:),'EdgeColor','c'); 
      if parameters.showNumbers
         text(box_coords(n,1),box_coords(n,2),num2str(n),'color','w');
      end
    end
end


   
