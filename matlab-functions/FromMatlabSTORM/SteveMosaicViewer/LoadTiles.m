function [mosaicImage,mosaicPars] = LoadTiles(varargin)
% [mosaicImage,mosaicPars] = LoadTiles() 
%  -- Loads the steve mosaic currently in global variable stvfile and plots
%  the results in figure(1) for the first 1000 frames.   
% [mosaicImage,mosaicPars] = LoadTiles(mosaicFile) 

%%
global stvfile steveMosaic

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'showall', 'boolean', true}; 
defaults(end+1,:) = {'showbox', 'boolean', false};  % plot box over seed position 
defaults(end+1,:) = {'shrk', 'positive', 1}; 
defaults(end+1,:) = {'position', 'array', [0,0]};  % seed position
defaults(end+1,:) = {'showNumbers','boolean',true};
defaults(end+1,:) = {'N', 'positive', 2000};  % max number of tiles around seed position 
defaults(end+1,:) = {'frameSize', 'positive', 1024}; 
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------

if length(varargin) == 0
    mosaicFile = stvfile;
    varInput = [];
else length(varargin) 
    mosaicFile = varargin{1};
    varInput = [];
end
if length(varargin) > 1
    varInput = varargin(2:end);
end

parameters = ParseVariableArguments(varInput, defaults, mfilename);

% default parameters
verbose = parameters.verbose; 
showall = parameters.showall;
position = parameters.position;
showbox = parameters.showbox;
shrk = parameters.shrk; 
N = parameters.N; 



%% Main Function

[mosaicFolder,mosaicName] = fileparts(mosaicFile);
mosaicFolder = [regexprep(mosaicFolder,'\','/'),'/']; 


% Read in .mat files and get list
%-----------------------------------
qstart = dir([mosaicFolder,'quickstart.mat']);

if isempty(qstart); 
    mtiles = dir([mosaicFolder,mosaicName,'*.mat']);
    M = length(mtiles); 
    if M==0
        disp(['found ',num2str(M),' tiles in folder ',mosaicFolder]); 
        error('no .mat mosaic files found');
    end
    if verbose
        disp(['found ',num2str(M),' tiles in folder ',mosaicFolder]);
    end

    xu = zeros(M,1);
    yu = zeros(M,1);
    xp = zeros(M,1);
    yp = zeros(M,1); 
    z = zeros(M,1); 
    mag = zeros(M,1);

    for m=1:M
         load([mosaicFolder,mtiles(m).name]);
             xp(m) = x_pix;
             xu(m) = x_um;
             yp(m) = y_pix;
             yu(m) = y_um;
             mag(m) = magnification;
             z(m) = zvalue;
             if verbose
                disp(['sorting data... ',num2str(100*m/M,3),'%']);
             end
    end

    save([mosaicFolder,'quickstart.mat'],'xu','yu','xp','yp','z','mag','mtiles','M'); 
else
    load([mosaicFolder,'quickstart.mat']); 
end
    
%%
N = min(N,M); 
[frames,~] = knnsearch([xu,yu],position(1,:),'k',N);

% get dimensions of grid to plot:
if showall
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
mosaicImage = zeros(ys,xs,1,'uint16');

for k = 1:length(frames); 
    i = frames(k);
    load([mosaicFolder,mtiles(i).name]);
    im = squeeze(data);

    % for plotting the 4x images, need to scale up:
    if showall
        if mag(i) ~= 1
            im = imresize(im,1/mag(i));
        end 
    end
    
 % rescale image if requested
     if shrk ~= 1
        im = imresize(im,1/shrk); 
     end
    [h,w,~] = size(im);
    h2 = round(h/2); % avoid integer operands for colon operator warning
    w2 = round(w/2);    
    y1 = Y(i)+1-h2;
    y2 = Y(i)+h2;
    x1 = X(i)+1-w2;
    x2 = X(i)+w2;
    
    
  try  
     mosaicImage(y1:y2,x1:x2) = im';
  catch
      mosaicImage(y1:y2-1,x1:x2-1) = im';
  end
  %   mosaicImage(x1:x2,y1:y2) = fliplr( flipud(im) );

    if verbose
        if rem(k,5) == 0
          disp(['rebuilding mosaic ', num2str(100*k/length(frames),3), '% done']);
        end
    end
%     figure(1); 
%     imagesc(imresize(mosaicImage,1)); 
%     pause(.3); colormap gray;
end

%--------- show plot
clf;
imOut = imadjust(mosaicImage,stretchlim(mosaicImage,0.001));
imagesc(imOut); 
colormap(gray(2^8));
set(gcf,'color','w');
axis image;
steveMosaic = imOut;
disp('saved mosaic image in global variable steveMosaic'); 

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

if showbox
    hold on;
    for n=1:numPos
      rectangle('Position',box_coords(n,:),'EdgeColor','c'); hold on;
      if parameters.showNumbers
        text(box_coords(n,1),box_coords(n,2),num2str(n),'color','w');
      end
    end
end


   
