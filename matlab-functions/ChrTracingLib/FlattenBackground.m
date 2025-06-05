function [imTiles,bkd] = FlattenBackground(imTiles,varargin)
% input is a cell array of 2D images 
%   OR a matrix of 2D images HxWxN where N is number of images
% 
% linear noise model. image = p*I + c, where p is the illumination profile
% and c is a constant pixel background. 

defaults = cell(0,3);
defaults(end+1,:) = {'smoothRes','positive',1}; % .1
defaults(end+1,:) = {'removeZeros','boolean',true}; % Translated images are padded with zeros, which can lead to bright borders.  
defaults(end+1,:) = {'bkdPerc','fraction',.01};
defaults(end+1,:) = {'bkd','freeType',[]}; % map of illumination background
defaults(end+1,:) = {'camBkd','freeType',[]}; %empty for zero, 'auto' for compute, otherwise should be a map  
defaults(end+1,:) = {'backgroundCorrect',{'file','median','medianEdge','none','removeData'},'median'};
defaults(end+1,:) = {'theta','positive',1.2}; % .1  % threshold for "removeData" mode    
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'strength','fraction',0.5};
pars = ParseVariableArguments(varargin,defaults,mfilename);

bkd = 0; % default value for assignment

if strcmp(pars.backgroundCorrect,'none')
    return
else    
    if iscell(imTiles)
        nTiles = length(imTiles);
        [h,w] = size(imTiles{1});
        imType = class(imTiles{1});
        imStack = double(cat(3,imTiles{:}));
    else
        imType = class(imTiles);
        imStack = double(imTiles);
        [h,w,nTiles] = size(imTiles);
    end
    if ischar(pars.camBkd)
       camBkd = min(imStack,[],3);
    elseif size(pars.camBkd) > 1
        camBkd = pars.camBkd;
    elseif isempty(pars.camBkd)
        camBkd = zeros(h,w);
    else
        cprintf([1 0 0],'did not recognize ');
        cprintf([1 0 0],pars.camBkd);
    end
    if pars.removeZeros
       low = quantile(imStack(:),pars.bkdPerc);
    end
    
% follow selected program
    switch pars.backgroundCorrect 
        case 'file'
            if isempty(pars.bkd)
                error('must supply a background image to use this option');
            else
                bkd = pars.bkd;
            end
        case 'median'
            if pars.removeZeros
                imStack(imStack==0) = low; 
            end
            imStack = imStack - repmat(camBkd,1,1,nTiles);
            bkd = nanmedian(imStack,3);    
        case 'removeData'
             [v,x] = hist(double(imStack(:)),0:20:2^15);
             % figure(2); clf; bar(x,v);
             [pks,loc] = findpeaks(v,x,'MinPeakDistance',10,'MinPeakHeight',1000,'MinPeakProminence',100);
             theta = pars.theta*loc(1);
             bkdStk = nan(h,w,nTiles);
            for f=1:nTiles
                im = imStack(:,:,f); % imTiles{f};
                if pars.removeZeros
                    im(im==0) = low; 
                end
                bw = false(size(im));
                bw(im >theta) = true;
%                 figure(3); clf; 
%                 subplot(1,2,1); imagesc(im); caxis([0,1500]);
%                 subplot(1,2,2); imagesc(bw); colorbar;
                im2 = double(im);
                im2(bw) = nan;
                bkdStk(:,:,f) = im2;
                pause(.1);
            end
            bkd = nanmedian(bkdStk,3);
            bkd = regionfill(bkd,isnan(bkd)); % interpolate missing region

        case 'medianEdge'
            % median but ignoring edge values
            topEdge = imStack(1:20,:,:);
            lftEdge = imStack(:,1:20,:);
            botEdge = imStack(end-19:end,:,:);
            rgtEdge = imStack(:,end-19:end,:);
            edgBkd = median( cat(1, topEdge(:),lftEdge(:),botEdge(:),rgtEdge(:)) );
            imStack2 = imStack;
            imStack2(imStack2<edgBkd) = nan;
             bkd = nanmedian(imStack2,3);   
             bkd(isnan(bkd)) = edgBkd;
             bkd = imgaussfilt(bkd,100);
             figure(10); clf; imagesc(bkd);
        otherwise
            disp(['correction option ', backgroundCorrect ,' not recognized.']);
    end



    
    if pars.smoothRes ~=1
        bkd = imresize(imresize(bkd,pars.smoothRes),[h,w]);
    end
    
    sc = quantile(bkd(:),pars.strength);
    if pars.showPlots
        imOut = cast(sc*imStack(:,:,1)./bkd,imType);
        figure(10); clf;
        subplot(1,3,1); imagesc(IncreaseContrast(bkd,'high',.999)); colorbar; title('background');
        subplot(1,3,2); imagesc(IncreaseContrast(imStack(:,:,1),'high',.999)); colorbar; title('original');
        subplot(1,3,3); imagesc(IncreaseContrast(imOut,'high',.999)); colorbar;  title('corrected');
        pause(1);
    end   
    % apply linear background subtraction
    if iscell(imTiles)
        for t=1:nTiles
            imTiles{t} = cast(sc*imStack(:,:,t)./bkd,imType);
        end
    else
        for t=1:nTiles
            imTiles(:,:,t) = cast(sc*imStack(:,:,t)./bkd,imType);
        end
    end
end


% imFixed = cast(quantile(bkd(:),.5)*double(imOrig)./bkd,imType);