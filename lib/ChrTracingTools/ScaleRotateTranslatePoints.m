function [ptsOut,imageTiles] = ScaleRotateTranslatePoints(ptsIn,alignValues,varargin)
% ptsIn - Nx2 vector of (x,y) points that indicate upper left coords of an
%        image.
% theta - scalar, rotation angle in degrees (not radians)
%       - struct - alignValues struct containing fields theta, xshift,
%                yshift, possibly rescale (ignored) 
% 'center' - default [0,0]. Compute rotation around this point
% 'shifts' - default [0,0]. Translate by this after the rotation
% 'imageTiles' - default empty. apply rotation translation to upper-left
%           tile position and image tile pairs (ul imageTile pairs).
% RotateTranslatePoints(uls,alignValues,'center',[],'imageTiles',imageTiles)
% 
% ScaleRotateShift, to simplify order of operations and allow more of them.
% 
defaults = cell(0,3);
defaults(end+1,:) = {'center','array',[0,0]};
defaults(end+1,:) = {'shifts','array',[0,0]};
defaults(end+1,:) = {'imSize','array',[0,0]};
defaults(end+1,:) = {'imageTiles','cell',{}};
defaults(end+1,:) = {'upperLeft','boolean',true};
defaults(end+1,:) = {'invert','boolean',false};
defaults(end+1,:) = {'opOrder',{'TR','RT'},'RT'}; % rotate, then translate
% should generalize to allow rescaling
pars = ParseVariableArguments(varargin,defaults,mfilename);

% unpack alignValues
scale1 = alignValues.rescale;
theta1 = alignValues.theta;
xshift1 = alignValues.xshift;
yshift1 = alignValues.yshift;
if isfield(alignValues,'rescale2')
    scale2 = alignValues.rescale2;
else
    scale2 = 0;
end
if isfield(alignValues,'theta2')
theta2 = alignValues.theta2;
else
    theta2 = 0;
end
if isfield(alignValues,'xshift2')
    xshift2 = alignValues.xshift2;
    yshift2 = alignValues.yshift2;
else
    xshift2 = 0;
    yshift2 = 0;
end

% don't accept multistep changes at present
if isfield(alignValues,'theta2')
    if alignValues.theta2 ~= 0
        error('does not currently do nested. Call function sequential');
    end
end

% check dimension of inputs
[nPts,dim] = size(ptsIn);
if dim~=2
    error('pts in must be nx2');
end

% 
imageTiles = pars.imageTiles;
if ~isempty(pars.imageTiles)
    [h_i,w_i,~] = size(imageTiles{1});
    pars.imSize = [h_i,w_i];
end
if sum(pars.imSize)~=0 
    cent1 = pars.center + [pars.imSize]/2*scale1;
    cent2 = pars.center + [pars.imSize]/2*scale1*scale2;
else
    error('please pass an image size or an image');
end

% scaling matrix
S1 = [scale1   0     0;
     0       scale1  0;
     0         0      1];

 % rotation matrix around 0,0
R = [cosd(-theta1) -sind(-theta1)  0
     sind(-theta1)  cosd(-theta1)  0
          0           0      1];

% Shift to center matrix (to rotate arond point 'cent').
P = [1 0 -cent1(1)
     0 1 -cent1(2)
     0 0 1];

R1 = P\R*P; 
 
 % tranlation matrix
T1 = [1 0 xshift1
     0 1 yshift1
     0 0 1     ];
 
 % scaling matrix
S2 = [scale2   0     0;
     0       scale2  0;
     0         0      1];

 % rotation matrix around 0,0
R = [cosd(-theta2) -sind(-theta2)  0
     sind(-theta2)  cosd(-theta2)  0
          0           0      1];

% Shift to center matrix (to rotate arond point 'cent').
P = [1 0 -cent2(1)
     0 1 -cent2(2)
     0 0 1];

R2 = P\R*P; 
 
 % tranlation matrix
T2 = [1 0 xshift2
     0 1 yshift2
     0 0 1     ];
 
 

data = [ptsIn'; ones(1,nPts)];
% first scale rotate shift
centMatrix = [repmat(cent1',[1,nPts]); ones(1,nPts)];
rescaled = centMatrix-scale1*(centMatrix  - data);
rotated = R1*rescaled;
translated = T1*rotated;

% second scale rotate shift
centMatrix = [repmat(cent2',[1,nPts]); ones(1,nPts)];
rescaled = centMatrix-scale2*(centMatrix  - translated);
rotated = R2*rescaled;
translated = T2*rotated;

        
% compute updated UL positions
ptsOut = translated(1:2,:)';

% apply scaling and rotation to image;
if ~isempty(imageTiles)
    [nTiles,nChns] = size(imageTiles);
    for n=1:nTiles
        for c=1:nChns
            im = imageTiles{n,c};
            im = ScaleRotateIm(im,alignValues);
            imageTiles{n,c} = im;
        end
    end
end
