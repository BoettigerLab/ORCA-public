function [fidiSpots,dataSpots,pars] = ChrTracer3p3_CropAllSpots(rawDataNames,regData,eTable,spots,varargin)
% 
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% 'parameters',pars
% A structure containing pars.fiducialFrames, and pars.dataFrames
% Both are 4D matrices, nRows x nCols x nStacks x nHybes 
% 
% -------------------------------------------------------------------------
% Optional inputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% fidiSpots{s}{h}
% 
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
% called by ChrTracer

% supress some unnecessary warnings. 
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off','MATLAB:prnRenderer:opengl');

defaults = cell(0,3);
% current spot 
defaults(end+1,:) = {'currentSpot','integer',1};
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'saveFolder','string',''};
% FOV parameters
defaults(end+1,:) = {'fov', 'integer', 0};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
% advanced / exploratory
defaults(end+1,:) = {'downSampleZ','integer',1}; % 
defaults(end+1,:) = {'checkOffset','boolean',false}; % 
defaults(end+1,:) = {'troubleshoot','boolean',false}; % 
pars = ParseVariableArguments(varargin, defaults, mfilename);




%%
tic;

% optional parameters
boxWidth = pars.boxWidth; % 16;
numHybes = size(rawDataNames,1);


%% parse info file 

% ------ Identify fiducial channel and parameter channel ------
h  = 1; % Moved outside loop, not a current working option to allow different number of data channels per hyb
[isFidChannel,frameChannels,~,~,dataChns] = GetFidChnFromTable(eTable,'hyb',h);
nBufFrames = eTable.bufferFrames(h);
totFrames  = eTable.totalFrames(h);
skipFrames = false(totFrames,1);
skipFrames([1:nBufFrames,totFrames-nBufFrames+1:totFrames]) = true;
isFidChannel(skipFrames) = false; % could move out of loop
numDataChns = length(dataChns);

isDatChannel = cell(numDataChns,1);
for n=1:numDataChns
    isCurrDat = StringFind(frameChannels,dataChns{n},'boolean',true); 
    isDatChannel{n} =  ~skipFrames & isCurrDat;
end
    

% =============  check zoffset
if pars.checkOffset    % check and correct off set to hybe 1
    zshift = ComputeHybeZoffset(rawDataNames,'parameters',pars);
end
%=====================
    
numSpots = size(spots,1);
% It is faster to load the data looping over hybes and getting all spots
% for each hybe. Later we will resort the data by looping over spots, which
% will be easier for analysis.
fidSpotsAllHybes = cell(numHybes,1);
datSpotsAllHybes = cell(numHybes,1);

for h=1:numHybes % parfor this  (doesn't help much, read-speed is limiting)
    dax = ReadDax(rawDataNames{h},'verbose',pars.veryverbose);
    dax = ApplyReg(dax,regData(h));
    [numRows,numCols,numZs] = size(dax); %#ok<ASGLU>
  
    fidSpotsInHyb = cell(numSpots,1);
    datSpotsInHyb = cell(numSpots,numDataChns);
    for s=1:numSpots        
        yi = round(max(1,spots(s,2)-boxWidth/2));
        ye = round(min(spots(s,2)+boxWidth/2,numRows));
        xi = round(max(1,spots(s,1)-boxWidth/2));
        xe = round(min(spots(s,1)+boxWidth/2,numCols));
        % multidimensional indexing
        fidMat = zeros(boxWidth+1,boxWidth+1,sum(isFidChannel));  % pad with zeros to ensure all data is same size
        fidMatIn = dax(yi:ye,xi:xe,isFidChannel);
        fidMat(1:size(fidMatIn,1),1:size(fidMatIn,2),1:size(fidMatIn,3)) = fidMatIn;
        fidSpotsInHyb{s} = fidMat;          
        % get data
        for n=1:numDataChns  
            datMat = zeros(boxWidth+1,boxWidth+1,sum(isDatChannel{n}));  % pad with zeros to ensure all data is same size
            datMatIn =  dax(yi:ye,xi:xe,isDatChannel{n});
            datMat(1:size(datMatIn,1),1:size(datMatIn,2),1:size(datMatIn,3)) = datMatIn;
            datSpotsInHyb{s,n} = datMat;
        end    
%         figure(2); 
%         subplot(1,2,1); imagesc(max(fidSpotsInHyb{s},[],3)); title('fid');
%         subplot(1,2,2); imagesc(max(datSpotsInHyb{s,n} ,[],3)); title('dat');
    end
    fidSpotsAllHybes{h}=fidSpotsInHyb;
    datSpotsAllHybes{h}=datSpotsInHyb;
end
tt = toc;
disp(['Loaded cropped data in ',num2str(tt/60),' min.']);
disp('Organizing hybs:');
tic

fidiSpots = cell(numSpots,1);
dataSpots = cell(numSpots,1);
for s=1:numSpots
    fidiSpots{s} = cell(numHybes,1);
    dataSpots{s} = cell(numHybes,numDataChns);
   for h=1:numHybes
       if pars.downSampleZ~=1
          im_fid = fidSpotsAllHybes{h}{s}(:,:,1:pars.downSampleZ:end); 
          im_dat = datSpotsAllHybes{h}{s,:}(:,:,1:pars.downSampleZ:end);
       else
          im_fid = fidSpotsAllHybes{h}(s); 
          im_dat = datSpotsAllHybes{h}(s,:);
       end
      fidiSpots{s}(h) = im_fid;
      dataSpots{s}(h,:) = im_dat; 
   end
end

%--------------
if pars.checkOffset
   for s=1:numSpots
       for h=1:numHybes
           imF = fidiSpots{s}{h};
           bkdFid = quantile(imF(:),.001);
           [h_i,w_i,z_i] = size(imF);
           if zshift(h,1) > 0 % shift hybe down towards data
               blnkF = bkdFid*ones(h_i,w_i,zshift(h,1),class(bkdFid));
               z_e = z_i-zshift(h,1); % +1?
               newImF = cat(3,blnkF,imF(:,:,1:z_e));
               newImD = cell(1,numDataChns);
               for n=1:numDataChns
                    imD = dataSpots{s}{h,n};
                    bkdDat = quantile(imD(:),.001);
                    blnkD = bkdDat*ones(h_i,w_i,zshift(h,1),class(bkdFid));
                    newImD{1,n} = cat(3,blnkD,imD(:,:,1:z_e));
               end
           else  % shift Hybe up
               blnkF = bkdFid*ones(h_i,w_i,-zshift(h,1),class(bkdFid));
               z_s = 1-zshift(h,1); % zshift is neg
               newImF = cat(3,imF(:,:,z_s:end),blnkF);
               newImD = cell(1,numDataChns);
               for n=1:numDataChns
                    imD = dataSpots{s}{h,n};
                    bkdDat = quantile(imD(:),.001);
                    blnkD =  bkdDat*ones(h_i,w_i,-zshift(h,1),class(bkdFid));
                    newImD{1,n} = cat(3,imD(:,:,z_s:end),blnkD);
               end 
           end
           fidiSpots{s}(h) = {newImF};
           dataSpots{s}(h,:) = newImD;
       end
   end
end
%---------
t=toc;
disp(['Finished organizing in ',num2str((t)/60),' min.']);
disp(['Finished Crop All in ',num2str((t+tt)/60),' min.']);


if pars.saveData 
    if ~isempty(pars.saveFolder)        
        fid4D = cell(numSpots,1);
        dat4D = cell(numSpots,1);
        for s=1:numSpots
            fid4D{s} = cat(4,fidiSpots{s}{:}) ;
            datTemp = dataSpots{s}';
            dat4D{s} = cat(4,datTemp{:}); % alternate chn1 and chn2
        end
        fid5D = cat(5,fid4D{:});
        dat5D = cat(5,dat4D{:});

        % mSize = cellfun(@size,fid4D,'UniformOutput',false);
        % mSize = cat(1,mSize{:})

        fullSaveName = [pars.saveFolder,'fov',num2str(pars.fov,'%03d'),'_FiduSpots.i5d'];
        WriteImage5D(fid5D,fullSaveName,...
            'fov',pars.fov,...
            'lociX',spots(:,1)',...
            'lociY',spots(:,2)');

        fullSaveName = [pars.saveFolder,'fov',num2str(pars.fov,'%03d'),'_DataSpots.i5d'];
        WriteImage5D(dat5D,fullSaveName,...
            'fov',pars.fov,...
            'lociX',spots(:,1)',...
            'lociY',spots(:,2)');
    else
        warning('specify a "saveFolder" in order to save data');
    end
end

%%

% 
% im4D = cat(4,fidiSpots{10}{:});
% figure(1); clf; PlotProjection4D(im4D,'projection','xz');
% figure(2); clf; PlotProjection4D(im4D,'projection','xy');

