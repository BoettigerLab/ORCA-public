function [spotCodeTable,maps,polys,spotCode] = DecodeSampleBarcodes(eTableXLS,analysisFolder,varargin)
%%
%  The experiment table should sepcify some of the hybes as dataType 'B'
%  for barcode.  LoadDaxFromEtable handles the data channel selection, the
%  detection and/or creation of max projections images, and drift
%  correction. 
% 
defaults = cell(0,3);
defaults(end+1,:) = {'scope','string','auto'};
defaults(end+1,:) = {'showFOVfig','integer',0}; 
defaults(end+1,:) = {'nHybs','integer',[]}; 
defaults(end+1,:) = {'nppXY','freeType',[]}; % will read from table eTableXLS info file by default
defaults(end+1,:) = {'overlayFig','integer',10}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'balanceBrightness','boolean',true}; 
defaults(end+1,:) = {'saveFigure','boolean',true}; 
defaults(end+1,:) = {'overwrite','boolean',true}; 
defaults(end+1,:) = {'showExtraPlots','boolean',false}; 
defaults(end+1,:) = {'recordQuality','boolean',false}; 
defaults(end+1,:) = {'codebook','string',''};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'maxProject','boolean',true}; % use full 3D stack or just max project.
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,'DecodeSampleBarcodes');
%%
% dataFolder = 'W:\2020-11-05_mES_dTAG_wo_barcodes_Hoxa\DNA_Expt\';

if pars.verbose
    disp('loading barcode images, please wait:');
end

[imMax,imInfo] = LoadDaxFromEtable(eTableXLS,...
    'hybNumber',inf,...
    'fov',inf,...
    'hybType','B',...
    'dataType','data',...
    'fixDrift',true,...
    'maxProject',pars.maxProject,...
    'driftFolder',analysisFolder,...
    'verbose',pars.verbose);

imMax0 = imMax; % backup for debugging only.

% remove empty spaces for skipped hybs
%     if barcodes are the last N hybs, there will be H blank cells which
%     were skipped for the readout hybs prior to the cellular barcodes. 
noData = cellfun(@isempty,imMax);
imMax(noData(:,1),:,:) = [];
imInfo(noData(:,1),:,:) = [];

% stack barcodes
[nH,nFOVs,nDatChns] = size(imMax);

% interleave data in the order it was aquired
%   this is the same convention used in ChrTracer3() and the experiment
%   table, and will make data matching easier
if nDatChns>1
   imMaxStk = cell(nH*nDatChns,nFOVs);
   k=0;
   for h=1:nH
       for c=1:nDatChns
           k=k+1;
          imMaxStk(k,:) = imMax(h,:,c);
       end
   end
   imMax = imMaxStk;
end

if size(imMax,1) <= 1
   warning('barcode folders must be designated as dataType "B" in the experiment table');
   error('Did not find multiple barcode folders, check your experiment table'); 
end

% determine scope settings from info file of first image
if isempty(pars.nppXY)
    scopePars = GetScopeSettings(imInfo{1,1,1});
    pars.nppXY = scopePars.nmPixXY;
end
%%

% remove data coded "barcode = 0", this symbol marks unused barcodes 
eTable = readtable(eTableXLS);
isB = strcmp(eTable.DataType,'B');
datPropTable = DataChnsFromTable(eTable(isB,:));
reads = datPropTable.readout;
idZero = find(reads==0);
imMax(idZero,:) = [];

%%
[nB,nFOVs] = size(imMax);
nHybs = pars.nHybs;
polys = cell(nFOVs,1);
maps = cell(nFOVs,1);
spotCodes = cell(nFOVs,1); % table - x,y,fov,barcode-value
try
for f=1:nFOVs
    if pars.showFOVfig
        figure(pars.showFOVfig); clf;
    end
    
     % ----load spot data in fov
%      spotXYtab = [];
%      spotTableFile = [saveFolder,'fov',num2str(f,'%03d'),'_selectSpots.csv'];
%      if exist(spotTableFile,'file')~=0
%         spotXYtab = readtable(spotTableFile);
%      end
%      spotXY = round([spotXYtab.locusX,spotXYtab.locusY]);
% 
%    --- load ORCA AllFits table     
     orcaTable = readtable([analysisFolder,'fov',num2str(f,'%03d'),'_AllFits.csv']);
     if isempty(nHybs)
         nHybs = max(orcaTable.readout(strcmp(orcaTable.dataType,'H')));
     end
     [polys{f},maps{f},spotData] = TableToPolymer(orcaTable,'bins',nHybs);
     spotXY = round(spotData/pars.nppXY); 
     nSpots = size(spotXY,1);
     spotCodes{f} = zeros(nSpots,3+nB);

    for b=1:nB        
         % record spot brightness
         [nRows,nCols] = size(imMax{b,f});
         for s=1:nSpots
             xi = spotXY(s,1);
             yi =  spotXY(s,2);
             xi = max([min([xi,nRows]),1]);
             yi = max([min([yi,nCols]),1]);
             if ~isnan(xi)
                 spotCodes{f}(s,1) = xi;
                 spotCodes{f}(s,2) = yi;
                 spotCodes{f}(s,3) = f; 
                 spotCodes{f}(s,3+b) = imMax{b,f}(yi,xi);
             end
         end 
         
         if pars.showFOVfig
             figure(pars.showFOVfig);
             subplot(1,nB,b); imagesc( imMax{b,f}); caxis([1e3,2^16]); % camera saturated
         end
    end
end

catch er
    warning(er.getReport);
    disp('debug here');
    
end

%%

spotCode = cat(1,spotCodes{:});
totSpots = size(spotCode,1);
spotGroup = zeros(totSpots,1);
contrast = nan(totSpots,1);

if ~isempty(pars.codebook)
    if ischar(pars.codebook)
        codebookTable = readtable(pars.codebook);
    else
        codebookTable = pars.codebook;
    end
    codebookBarcodes = codebookTable{1,2:end};
    codebookNames = codebookTable{2:end,1};
    codebook = codebookTable{2:end,2:end};
    
    
    
     codeScore = double(spotCode(:,4:end)); 
     codebook2 = codebook;
     % % OLD order of barcodes in spotCode is H1_750, H2_750 ... H1_647, H2_647     
     % chns = cellfun(@str2double,datPropTable.chn);
     % [chns,chnSort] = sort(chns,'descend');  % Not sure this sort step is robust to all variations of using barcodes
     % reads = datPropTable.readout(chnSort);
%      codebook2 = zeros(size(codebook,1),length(reads));
%      for b=1:length(codebookBarcodes)
%          cs = find( reads == codebookBarcodes(b));
%          for c=Row(cs)
%             codebook2(:,c) = codebook(:,b);
%          end
%      end
     
     % % just for troubleshooting, compare tables
%      codeTable2 = cat(1,reads',codebook2)
%      codeTable1 = codebookTable{:,2:end}
     
     % we want to balance the barcodes since not all are equally bright
     % Also the 750 images are systematically much dimmer than the 647
     % We'll try rescaling by 90% quantile.
     % All the data gets to saturate the camera somewhere, so using too
     % high a percentile doesn't achieve much balancing
     % Some barcodes may be very low abundance (though this should be
     % avoidable), so we don't want to rescale to the mean/median -- those
     % might all be legit blanks. Though this should be reduced in
     % intentional mixing experiments and with codebooks that use
     % individual barcodes mutliple times.  
     
     if pars.balanceBrightness
         mB=quantile(codeScore,.9); %
         % mB(reads==0) = inf; % removing the intentionally blank.  % Blanks are now handled above  
         codeScore = codeScore ./ repmat(mB,totSpots,1); 
         codeScore(codeScore>1) = 1;
     end
     % Now we 
     % figure(1); clf; imagesc(codeScore(randperm(totSpots,50),:));
     groupName = cell(totSpots,1);
     for s=1:totSpots
%          s = randperm(totSpots,1);
%           codeScore(s,:)
        [idx,dis] = knnsearch(codebook2,codeScore(s,:));
        spotGroup(s) = idx;
        contrast(s) = 1./dis;
        groupName{s} = codebookNames{idx};
     end
     %  nB = max(spotGroup);  %                              problematic  -- blanks are now handled above 

     % Troubleshooting
     %   check that all cell types were detected, plot abundances 
%         codebookNames
%         [uNames,ui] = unique(groupName)
%         grpCnt = hist(spotGroup,1:max(spotGroup));
%         figure(3); clf; barh(grpCnt);
%         set(gca,'YTickLabels',codebookNames)
%         xlabel('spot count'); set(gcf,'color','w');
%         % sort the rows of the codebook by abundance 
%         %  just make sure there are no uncanny patterns
%         [~,cntSort] = sort(grpCnt)
%         codebook2(cntSort,:)
%         % as expected, all the singles are the most abundant

     
    % convert to table
    x = spotCode(:,1);
    y = spotCode(:,2);
    fov = spotCode(:,3);
    spotCodeTable = table(x,y,fov,spotGroup,contrast,groupName);
     
else % OBSOLETE? 03/24/21


    % classify with max signal
    for s=1:totSpots
        [vs,is] = sort(spotCode(s,4:end),'descend'); 
        spotGroup(s) = is(1);
        contrast(s) = vs(1)./vs(2);
    end


    % figure(2); clf; hist(contrast,100);
    % figure(3); clf; hist(spotGroup,1:nB);

    %%
    if pars.balanceBrightness
        % Purpose
        % Some barcodes strip off poorly, such that the brightest cells
        % even after strand-displacement are on par with the dim cells
        % (e.g. spots on the edge of cells, spots on the edge of the FOV).
        % This is more prevelant with less abundant cell types (suggestion
        % adapter/dye is subsaturating). 
        % These bright hybs end up with bimodal distributions of brightness
        % combining the real stuff in one peak and the misclassified in
        % another. The aim of this balancing is to fix this
        % 
        % Method
        % Groups are initially assigned on brightest hyb
        % Then compute median brightness in each initial group assignment
        % (assumes that the median of a group is still the 'on target' set)
        % These barcode-strength data is then used to rescale the data
        % 
        % Notes: I'm not sure this is any more robust than top decile.
        if pars.showExtraPlots
            figure(2); clf;
            for b=1:nB
                figure(2); subplot(nB,2,2*(b-1)+1); 
                hist(contrast(spotGroup==b),[0:.5:10]);
                title(['contrast, median:' num2str(nanmedian(contrast(spotGroup==b)))]);
                xlim([0,10]);
                subplot(nB,2,2*b);
                hist(spotCode(spotGroup==b,3+b),0:2000:2^16); xlim([0,2^16]);
                title(['brightness, median: ',num2str(nanmedian( spotCode(spotGroup==b,3+b)))] );
            end
        end

        mB = zeros(1,nB);
        for b=1:nB
            mB(b) = nanmedian(spotCode(spotGroup==b,3+b));
        end

        spotCode(:,4:end) = spotCode(:,4:end) ./ repmat(mB,totSpots,1); 
        totSpots = size(spotCode,1);
        spotGroup = zeros(totSpots,1);
        contrast = nan(totSpots,1);
        % classify with max signal
        for s=1:totSpots
            [vs,is] = sort(spotCode(s,4:end),'descend'); 
            spotGroup(s) = is(1);
            contrast(s) = vs(1)./vs(2);
        end
        % convert to table
        x = spotCode(:,1);
        y = spotCode(:,2);
        fov = spotCode(:,3);
        spotCodeTable = table(x,y,fov,spotGroup,contrast);


        if pars.showExtraPlots
            figure(3); clf;
            for b=1:nB
                figure(3); subplot(nB,2,2*(b-1)+1); 
                hist(contrast(spotGroup==b),[0:.5:10]);
                title(['contrast, median:' num2str(nanmedian(contrast(spotGroup==b)))]);
                xlim([0,10]);
                subplot(nB,2,2*b);
                hist(spotCode(spotGroup==b,3+b),(0:2000:2^16)/2^16); xlim([0,1]);
                title(['brightness, median: ',num2str(nanmedian( spotCode(spotGroup==b,3+b)))] );
            end
        end

    end
end
%% compress these into matrices
% the fov's can still be separated with the spotCodes table
try
maps = cat(3,maps{:});
polys = cat(3,polys{:});
catch
   maps = StackMapCellArray(maps);
   polys = StackMapCellArray(polys);
end

%% 
isGood = true(1,nFOVs);
if pars.overlayFig > 0 
    if pars.saveFigure
       SetFigureSavePath(analysisFolder); 
    end
    figF = figure(pars.overlayFig); clf;
    subF = gca;
    for fov = 1:nFOVs
        
        if pars.saveFigure
           figName = [analysisFolder,'fov',num2str(fov,'%03d'),'_Demultiplex.png'];
           if exist(figName,'file')
              continue 
           end
        end
        
        subF.NextPlot = 'replace';
        cMap = GetColorMap('hsvCut',nB);
        ncImage = cat(3, imMax{:,fov});
        for b=1:nB
            ncImage(:,:,b) =uint16( 2^15*double(ncImage(:,:,b))./double(mB(b)) );
        end
        Ncolor(ncImage,'colormap',cMap);
        colormap(cMap);
        colorbar;
        isIn = spotCode(:,3)==fov;
        xy = spotCode(isIn,1:2);
        hold on; plot(xy(:,1),xy(:,2),'k.','MarkerSize',20)
        for b=1:nB  
            isIn = spotCode(:,3)==fov & spotGroup==b;
            isLow = isIn & contrast<1.5;
            xy = spotCode(isIn,1:2);
            xyL = spotCode(isLow,1:2);
            subF.NextPlot = 'add';
            plot(subF,xy(:,1),xy(:,2),'.','color',cMap(b,:),'MarkerSize',16);
            plot(subF,xyL(:,1),xyL(:,2),'+','color',cMap(b,:),'MarkerSize',16);
            
        end
        if pars.saveFigure
           figName = ['fov',num2str(fov,'%03d'),'_Demultiplex'];
           SaveFigure(figF,'name',figName,'formats',{'png','fig'},'overwrite',pars.overwrite); 
        end
        if pars.recordQuality
          answer = questdlg('Record quality:', ... % question
                            'Reply ', ...  % pop-up label
                              'Good','Bad','Exit','Good'); % op1 op2 op3 default
          if strcmp(answer,'Good')
              isGood(fov) = true;
          elseif strcmp(answer,'Bad')
              isGood(fov) = false;
          elseif strcmp(answer,'Exit')
              break
          end
        end
    end
    if pars.recordQuality
        disp('You recorded the following FOV as "Bad":');
        disp(find(~isGood));
    end
end