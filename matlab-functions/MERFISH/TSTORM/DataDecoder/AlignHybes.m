function [tform,warpUncert] = AlignHybes(fedPos,varargin)
% Aligns hybe images based on starting positions of beads.  Returns
% structure of tforms to warp each hybe. 
%
% tform = AlignHybes(fedPos,'maxD',30,'warp2Hybe1',false,'matchFig',[])
% 
%--------------------------------------------------------------------------
% Required parameters
% fedPos -- cell array of nx2 doubles listing the x and y positions of each
% of the n feducials found in each hybe.  
% 
%--------------------------------------------------------------------------
% Optional Parameters
% 
% maxD = 30;
% warp2Hybe1 = false;
% matchFig = [];
% 
%--------------------------------------------------------------------------

troubleshoot = false; 

maxD = 30;
warp2Hybe1 = false;
matchFig = [];
alignFig = [];
daxFile = [];
showPlots = true;
rawData = {};
%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 2
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'maxD'
                maxD = CheckParameter(parameterValue,'positive','maxD');
            case 'warp2Hybe1'
                 warp2Hybe1 = CheckParameter(parameterValue,'boolean','warp2Hybe1');
            case 'matchFig'
                matchFig = CheckParameter(parameterValue,'handle','matchFig'); 
            case 'alignFig'
                alignFig  = CheckParameter(parameterValue,'handle','alignFig'); 
            case 'daxFile'
                daxFile = CheckParameter(parameterValue,'array','daxFile');
            case 'rawData'
                rawData = CheckParameter(parameterValue,'cell','rawData');
            case 'showPlots'
                showPlots = CheckParameter(parameterValue,'boolean','showPlots');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%%

numHybes = length(fedPos);
rawBeadIm = zeros(256,256,numHybes,'uint16'); 
clrmap = jet(numHybes); 
tformNN1 = cell(numHybes,1);
tform = cell(numHybes,1);
warpedStarts = cell(numHybes,1);
warpUncert = cell(numHybes,1); 
warpErrors = NaN*ones(numHybes,1);
xl = cell(numHybes,1);
yl = cell(numHybes,1); 
r = 5; % radius for color wheels of aligned hybes

if showPlots && ~isempty(daxFile)
    try
    daxImage = ReadDax(daxFile,'endFrame',5,'verbose',false);
    catch er
        disp(er.getReport);
        daxImage = zeros(256,256,5,'uint16');
    end
    if isempty(alignFig)
       alignFig = figure; clf;
    else
        figure(alignFig);
    end
    imagesc(mean(daxImage,3)); colormap gray;
    hold on;
end

if isempty(matchFig) && showPlots
   matchFig = figure; clf; 
end
if isempty(alignFig)  && showPlots
    alignFig = figure; clf;
end

for i=1:numHybes
    if i>1
        hybe1 = fedPos{i-1};    
    else
        hybe1 = fedPos{1};
    end
    hybe2 = fedPos{i};
    
    if ~isempty(rawData)  && showPlots 
        try
        rawBeadImage(:,:,i) = ReadDax(regexprep(rawData{i},'_list.bin','.dax'),'endFrame',1,'verbose',false);
        figure(10); subplot(round(numHybes/4),4,i); imagesc(5*rawBeadImage(:,:,i));  colormap gray;
        % figure(3); clf; Ncolor(5*rawBeadImage,jet(numHybes));
        catch
        end
    end
    
    try
    if warp2Hybe1 % warp all beads to first hybe
        hybe1 = fedPos{1};
        [tform{i}, warpErrors] = Warp2BestPair(hybe1,hybe2,...
            'maxD',maxD,'fighandle',matchFig,'showPlots',showPlots); 
    else  % warp to previous hybe
        [tformNN1{i}, warpErrors] = Warp2BestPair(hybe1,hybe2,...
            'maxD',maxD,'fighandle',matchFig,'showPlots',showPlots);  
        tform{i} = maketform('composite',tformNN1{1:i});
    end
    catch er
        disp(er.getReport);
        warning(['failed to warp data for hybe ',num2str(i),'. Using previous positions']);
        if i==1
            tform{i} = maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
            tformNN1{i}= maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
        else
            tform{i} = tform{i-1}; % don't warp if you can't warp
            tformNN1{i}= maketform('affine',[1 0 0; 0 1 0; 0 0 1]); % don't warp if you can't warp
        end
    end
    
    numFeducials = length(warpErrors);
    disp(160*warpErrors(1:min(numFeducials,5))');
    warpUncert{i} = 160*warpErrors(1:min(numFeducials,5));
    
    [xw,yw] = tforminv(tform{i},hybe2(:,1),hybe2(:,2));
    warpedStarts{i} = [xw,yw];
    if showPlots
        figure(alignFig);
        plot(xw,yw,'.','color',clrmap(i,:)); hold on;  
        theta = pi/numHybes*i;
        xl{i} = reshape([xw-r*cos(theta),xw+r*cos(theta),NaN*ones(length(xw),1)]',length(xw)*3,1);
        yl{i} = reshape([yw-r*sin(theta),yw+r*sin(theta),NaN*ones(length(yw),1)]',length(xw)*3,1);
        
        
    end
end

if showPlots
    legend(num2str([1:numHybes]'),'Location','Best')
    set(gcf,'color','w'); 
    for i=1:numHybes
        plot(xl{i},yl{i},'color',clrmap(i,:),'LineWidth',2); hold on;  
    end
    xlim([0,256]); ylim([0,256]);
    xlabel('pixels'); ylabel('pixels');
end
