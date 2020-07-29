function [isFidChannel,frameChannels,channels,midFrame,dataChns] = GetFidChnFromTable(eTable,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'hyb', 'integer', 1};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin, defaults, mfilename);


    h = pars.hyb;
    channels = strsplit(regexprep(eTable.channels{h},'[^A-Za-z0-9,]',''),',');
    nBufFrames = eTable.bufferFrames(h);
    totFrames  = eTable.totalFrames(h);
    frameChannels = cell(totFrames,1);
    numChns = length(channels);
    fidChannel = strcmp(channels,num2str(eTable.fiducialChannel(h)));
    dataChns = channels(~fidChannel);
    fidChn = channels(fidChannel);
    % populate alternating channels
    for c=1:numChns
       nFr = length(frameChannels(c:numChns:end));
       frameChannels(c:numChns:end) = repmat(channels(c),nFr,1);
    end
   
    % if first channel is not specified as active or inactive, we assume it
    % is the old style inactive format. Thus we need to skip it. 
    if ~any(contains(eTable.Properties.VariableNames,'FirstFrameActive'))
        frameChannels = cat(1,'blnk',frameChannels); % first channel is blank and is a in the same color as the second channel. thereafter they alternate.   
        frameChannels = frameChannels(1:totFrames); % drop the last
    else
        if ~eTable.FirstFrameActive(1) % if first channel is specied as inactive skip it
            frameChannels = cat(1,'blnk',frameChannels); % add a blnk
            frameChannels = frameChannels(1:totFrames); % drop the last
        end
    end
    % 
    for f=1:nBufFrames
        frameChannels{f} = [frameChannels{f},'-m'];
        frameChannels{end+1-f} = [frameChannels{end+1-f},'-m'];
    end
    isFidChannel = StringFind(frameChannels,fidChn{1},'boolean',true);
    nFrames = totFrames -2*nBufFrames;
    midFrame = round(nFrames/numChns);