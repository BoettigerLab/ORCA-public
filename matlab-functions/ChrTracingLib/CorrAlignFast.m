function [alignValues,pars] = CorrAlignFast(im1,im2,varargin)
%  [alignValues,parameters] = CorrAlignFast(Im1,Im2)
% Inputs
% {'maxSize', 'positive', 200}; image size to start acceleration
% {'showplot', 'boolean', false}; show image of before and after  
% 
% Compute xshift and yshift to align two images based on maximizing
% cross-correlation.  
% 
% Updates 2019-10-23
% unless we are only doing translation, it is not possible to collapse the
% coarse and fine alignment steps.  
% We could in future write the combination function to reduce 
% scale, rotate, shift, scale, rotate, shift, to a single set of
% transactions (it will need to know the rotation center to use for the
% fine align). Meantime, lets rewrite the processing functions to do the 2x
% computation. For clarity and backwards compatibility (though not 
% elegance), chose to explicit represent this as 2 step process with a
% xshift2, theta2, etc.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment. set to 'inf' to skip to step alignment
defaults(end+1,:) = {'fineBox', 'freeType', []};  % perform fine scale alignment using a box of this size around the brightest point.
defaults(end+1,:) = {'fineUpsample', 'positive', 1};  % 
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % -10:1:10
defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 30};
defaults(end+1,:) = {'fineAngles','float',0}; % -1:.1:1
defaults(end+1,:) = {'fineScales','float',1}; % 0.95:0.01:.1.05
defaults(end+1,:) = {'fineCenter','array',[0,0]}; % x0 y0, [0,0] = use brightest
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'fastDisplay', 'boolean', true};
defaults(end+1,:) = {'displayWidth', 'integer', 500};
defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
defaults(end+1,:) = {'minFineImprovement', 'float', .1}; % the correlation coefficient of the fine scale alignment must be at least this fold greater than that of the coarse scale alignment to be used. Otherwise only coarse alignment is used
defaults(end+1,:) = {'label1', 'string', ''};
defaults(end+1,:) = {'label2', 'string', ''};
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'two 2D image matrices are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
pars = ParseVariableArguments(varargin, defaults, mfilename);

% %%
% pars = ParseVariableArguments([], defaults, mfilename);
% im1 = IncreaseContrast(imageTiles{1},'low',.7,'high',.999); figure(10); clf; imagesc(im1);
% im2 = TranslateImage(imrotate(im1,0,'bilinear','crop'),-120,300); figure(2); clf; imagesc(im2);
%%
if isempty(pars.fineBox)
    pars.fineBox = round(pars.maxSize/2);
end
[H,W] = size(im1);
relSpeed = sqrt(H*W)/pars.maxSize;
if relSpeed > 30 && pars.verbose
    warning(['downsampling ',num2str(round(relSpeed)),...
        ' fold, results may be inaccurate.  Try increasing maxSize to improve accuracy. Current maxSize=',num2str(pars.maxSize)]);
end
if relSpeed > 1
    % downsample and compute coarse correction
    % im1s = imresize(im1,1/relSpeed);
    % im2s = imresize(im2,1/relSpeed);   
       S = [1/relSpeed   0   0;
             0   1/relSpeed  0;
             0       0       1];
      im1s = imwarp(im1,affine2d(S) ,'OutputView',imref2d(  round(size(im1)/relSpeed))  );
      im2s = imwarp(im2,affine2d(S) ,'OutputView',imref2d(  round(size(im2)/relSpeed))  );
    
    if pars.showExtraPlot
        figure(31); clf;
    end
    if length(pars.scales) > 1 && sum(pars.scales~=1)
        % scales = 1-(1-pars.scales)/relSpeed
        scales = pars.scales;
    else
        scales = 1;
    end
    [alignValues1,car1]=CorrAlignRotateScale(im1s,im2s,...
        'parameters',pars,...
        'maxShift',pars.maxShift/relSpeed,...
        'showplot',pars.showExtraPlot,...
        'scales',scales,...
        'verbose',false);
    im2b = ScaleRotateShift(im2,alignValues1,'rescaleShifts',relSpeed);
    % Fine scale correction using a subset of the (now coarse aligned) image 
    if sum(pars.fineAngles) ~= 0
        % Fine scale rotation requires using the center field of view
        cx = round(H/2);
        cy = round(Y/2);
    elseif sum(pars.fineCenter) == 0
        % No rotation requested, we can pick a region with lots of data:
        % Downsampling also finds a region with substantial
        % data and not just a single bright pixel.
        subpix = 5; % 
        speedScale = max([0.01,1/(2*relSpeed*subpix)]);
        im3 = double(imresize(im1,speedScale)).*double(imresize(im2b,speedScale)); %  
%         figure(1); clf; subplot(1,3,1); imagesc(imresize(im1,speedScale)); colorbar;
%         subplot(1,3,2); imagesc(imresize(im2b,speedScale)); colorbar;
%         subplot(1,3,3); imagesc(im3); title('comb'); colorbar;
        [~,indmax] =  max(im3(:));
        [cy,cx] = ind2sub(size(im3),indmax);
        cy = round(cy/speedScale);  cx = round(cx/speedScale);
      
    else
        cx = pars.fineCenter(1); 
        cy = pars.fineCenter(2); 
    end
    x1 = max(1,cx-pars.fineBox);
    x2 = min(W,cx+pars.fineBox);
    y1 = max(1,cy-pars.fineBox);
    y2 = min(H,cy+pars.fineBox);
    pars.fineBoxCoords = table(x1,x2,y1,y2); 
    im1z = im1(y1:y2,x1:x2);
    im2z = im2b(y1:y2,x1:x2);
   %  figure(1); clf; imagesc(im3);
   %  figure(1); clf; imagesc(im1); hold on; rectangle('position',[x1,y1,x2-x1,y2-y1],'EdgeColor','w');
    if pars.showExtraPlot
        figure(32); clf;
    end
    [alignValues2,car2]=CorrAlignRotateScale(im1z,im2z,...
        'parameters',pars,...
        'angles',pars.fineAngles,...
        'scales',pars.fineScales,...
        'maxShift',round(relSpeed)+pars.fineMaxShift,...
        'showplot',pars.showExtraPlot,...
        'upsample',pars.fineUpsample);
    if car2.corrPeak < pars.minFineImprovement*car1.corrPeak
        alignValues2.xshift = 0; 
        alignValues2.yshift = 0;
        alignValues2.theta = 0;
        fineAlignFail = true;
        if pars.verbose
           disp('high res alignment failed, using only coarse alignment'); 
        end
    else
        fineAlignFail = false;
    end
    % save results for export
    % note: the order of events matter, so we can't just add the shifts and
    % the thetas.
    
    alignValues.xshift = round(relSpeed*alignValues1.xshift);
    alignValues.xshift2 = alignValues2.xshift;
    alignValues.yshift = round(relSpeed*alignValues1.yshift);
    alignValues.yshift2 = alignValues2.yshift;
    alignValues.theta = alignValues1.theta;
    alignValues.theta2 = alignValues2.theta;
    alignValues.rescale = alignValues1.rescale; % 1-relSpeed*(1-alignValues1.rescale); 
    alignValues.rescale2 = alignValues2.rescale;
    pars.corrPeak = car1.corrPeak;
    pars.corrPeak2 = car2.corrPeak;
    % plot results
    if pars.showplot
        if pars.showExtraPlot
            figure(33); clf;
        end
        if pars.fastDisplay
           displaySize = [pars.displayWidth*H/W,pars.displayWidth];
           im1_show = imresize(im1,displaySize);
           im2_show = imresize(im2,displaySize);
           im2b_show = imresize(im2b,displaySize);
        else
           im1_show = im1;
           im2_show = im2;
           im2b_show = im2b;
        end
        im2c_show = ScaleRotateShift(im2,alignValues); % will zoom in
        im2c_show = im2c_show(y1:y2,x1:x2);

        
        if ~isempty(pars.label1)
            labelIm = [pars.label1,'-red ',pars.label2,'-cyan'];
        else
            labelIm = '';
        end
        
        im3 = cat(3,im1_show,im2_show);
        sp1 = subplot(2,2,1); 
        im_sp1 = Ncolor(IncreaseContrast(im3,'high',.999));
        imagesc(sp1,1:W,1:H,im_sp1);
        title([labelIm,'  original']);
        im3 = cat(3,im1_show,im2b_show);
        sp3 = subplot(2,2,3);
        im_sp3 = Ncolor(IncreaseContrast(im3,'high',.999));
        imagesc(sp3,1:W,1:H,im_sp3);
        title({...
            labelIm,...
            ['coarse correction, ',num2str(relSpeed,2),' pixel accur.'],...
            ['xshift=',num2str(round(relSpeed*alignValues1.xshift)),...
            ' yshift=',num2str(round(relSpeed*alignValues1.yshift)),...
            ' rot=',num2str(alignValues1.theta),...
            ' rescale=',num2str(alignValues1.rescale),... 1-relSpeed*(1-alignValues1.rescale)
            ' peak=',num2str(car1.corrPeak)]});
        
        % im3 = cat(3,im1_show,im2c);
        im3 = cat(3,im1(y1:y2,x1:x2),im2c_show);
        sp4 = subplot(2,2,4); 
        im_sp4 = Ncolor(IncreaseContrast(im3,'high',.999));
        imagesc(sp4,(x1:x2),(y1:y2),im_sp4);
        % imagesc(sp4,relSpeed*(x1:x2),relSpeed*(y1:y2),im_sp4);
        % im3 = cat(3,im1,im2b); % for troubleshooting. 
       %  subplot(2,2,2); Ncolor(IncreaseContrast(im3,'high',.999));
        if ~fineAlignFail
            title({...
                labelIm,...
                ['fine correction, ',num2str(1/pars.fineUpsample,2),' pixel accur.'],...
                ['add. xshift=',num2str(alignValues.xshift2),...
                ' yshift=',num2str(alignValues.yshift2),...
                ' rot=',num2str(alignValues.theta2),...
                ' rescale=',num2str(alignValues.rescale2),...
                ' peak=',num2str(car2.corrPeak)]});
        else
           title({...
                labelIm,...
                ['FAILED fine correction, ',num2str(1/pars.fineUpsample,2),' pixel accur.'],...
                ['add. xshift=',num2str(alignValues.xshift2),...
                ' yshift=',num2str(alignValues.yshift2),...
                ' rot=',num2str(alignValues.theta2),...
                ' rescale=',num2str(alignValues.rescale2),...
                ' fine peak=',num2str(car2.corrPeak)]});
        end
    end
else
    [alignValues,pars] = CorrAlignRotateScale(im1,im2,'parameters',pars);
    alignValues.xshift2 = 0;
    alignValues.yshift2 = 0;
    alignValues.theta2 = 0;
    alignValues.rescale2 = 0;
end


