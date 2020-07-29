function dax = CroppedSpotsToDax(spotsCell,varargin)
% dax = CroppedSpotsToDax(spotsCell,varargin)
% defaults
% defaults(end+1,:) = {'zStart', 'integer', 10/2+1};
% defaults(end+1,:) = {'zStop', 'integer', 10/2};
% defaults(end+1,:) = {'dim','positive',3};
% defaults(end+1,:) = {'align','boolean',false}; 
% defaults(end+1,:) = {'goodHybes','array',[]}; 
% -------------------------------------------------------------------------
% Alistair Boettiger
% boettiger@stanford.edu
% Feb 10th, 2017
% Copyright CC BY NC
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'zStart', 'integer', 10/2+1};
defaults(end+1,:) = {'zStop', 'integer', 10/2};
defaults(end+1,:) = {'dim','positive',3};
defaults(end+1,:) = {'align','boolean',false}; 
defaults(end+1,:) = {'goodHybes','array',[]}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'increaseContrast','boolean',true}; 
defaults(end+1,:) = {'method',{'maxProject','zStack'},'maxProject'}; 
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% parse parameters 
ss = parameters.zStart ;% start Zscan (after discarded frames).  
es = parameters.zStop; % end Zscan 
numSpots = length(spotsCell);
if isempty(parameters.goodHybes)
    parameters.goodHybes = ones(numSpots,1);
end
[v,w,z,numHybes] = size(spotsCell{1});



if nargout>1
    alignSpotCell = spotCell;
end

% -------------------------------------------------------------------------
% Main Loop
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Method 1: max project
% -------------------------------------------------------------------------
if parameters.dim == 3
    dax = zeros(v+1,w+1,(numHybes+1)*numSpots,'uint16');
else
    zFrames = length(ss:z-es);
    dax = zeros(w+1,zFrames,(numHybes+1)*numSpots,'uint16');
end

switch parameters.method
    case 'maxProject'
        for s=1:numSpots
            if parameters.verbose
                disp(['processing spot ',num2str(s),' of ',num2str(numSpots)]);
            end
            im = squeeze(max(spotsCell{s}(:,:,ss:end-es,:),[],parameters.dim));
            im = permute(im,[2,1,3]); % standardizing who's x and y
            if parameters.increaseContrast
                im = IncreaseContrast(im,'low',0,'high',1)/uint16(4); % avoid saturating (at least the viewer saturates)
            end
            imOut = im; % just add the VxWxHybes block into the dax
            if parameters.align  
                % loop through h Hybes and coor align them to the previous good hybe  
                for h=2:numHybes
                    alignHybe = h-1; 
                    if ~parameters.goodHybes(h) % if current hybe is bad, skip it
                        continue
                    end
                    good = parameters.goodHybes(alignHybe); % if previous hybe is blank/bad, we'll use the one before it 
                    while ~good && alignHybe > 0
                        alignHybe = alignHybe-1;
                        good = parameters.goodHybes(alignHybe);
                    end
                    h1 = imOut(:,:,alignHybe);
                    h2 = im(:,:,h);
                    [xshift,yshift] = CorrAlign(h1,h2,'showplot',false,'upsample',4);
                    h2s = TranslateImage(h2,xshift,yshift,'upsample',4);
                    imOut(:,:,h) = h2s;
                    if nargout > 1
                        alignSpotCell{s}(:,:,:,h) = TranslateImage(alignSpotCell{s}(:,:,:,h),xshift,yshift,'upsample',4);
                    end
                end
            end
            [v,w,~] = size(imOut); % not all spots have same dimensions if they were near the image edge
            dax(1:v,1:w,(s-1)*(numHybes+1)+1:s*(numHybes+1)-1) = imOut;  % skips a frame between each spot 
        end

    % -------------------------------------------------------------------------
    % Method 2: whole z-stack
    % -------------------------------------------------------------------------    
    case 'zStack'
        dax = zeros(v+1,w+1,numHybes*numSpots*z,'uint16');
        for s=1:numSpots
            [v,w,z,numHybes] = size(spotsCell{s});
            spotsCell{s}(:,:,:,~parameters.goodHybes) = spotsCell{s}(:,:,:,~parameters.goodHybes)*uint16(0); % set bad hybes to zero
            dax(1:v,1:w,numHybes*z*(s-1)+1:numHybes*z*s) = reshape(spotsCell{s},[v,w,z*numHybes]);
        end
end
   



%%

%  % examples
% daxCy3_xy = CroppedSpotsToDax(cy3Spots,'zStart',ss,'zStop',es,'dim',3);
% daxName_xy = ['roi',num2str(d,'%03d'),'_spot001-',num2str(length(spots),'%03d'),'_cy3_xy'];
% WriteDax(daxCy3_xy,'folder',saveFolder,'daxName',daxName_xy);
% 
% daxCy3_xz = CroppedSpotsToDax(cy3Spots,'zStart',ss,'zStop',es,'dim',2);
% daxName_xz = ['roi',num2str(d,'%03d'),'_spot001-',num2str(length(spots),'%03d'),'_cy3_xz'];
% WriteDax(daxCy3_xz,'folder',saveFolder,'daxName',daxName_xz);
% 
% daxCy5_xy = CroppedSpotsToDax(cy5Spots,'zStart',ss,'zStop',es,'dim',3);
% daxName_xy = ['roi',num2str(d,'%03d'),'_spot001-',num2str(length(spots),'%03d'),'_cy5_xy'];
% WriteDax(daxCy5_xy,'folder',saveFolder,'daxName',daxName_xy);
% 
% daxCy5_xz = CroppedSpotsToDax(cy5Spots,'zStart',ss,'zStop',es,'dim',2);
% daxName_xz = ['roi',num2str(d,'%03d'),'_spot001-',num2str(length(spots),'%03d'),'_cy5_xz'];
% WriteDax(daxCy5_xz,'folder',saveFolder,'daxName',daxName_xz);