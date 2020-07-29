function [distMap,polymer,spotProps] = TableToDistMap(spotTable,varargin)
% distMap = TableToDistMap(bxcTable)
% sorts by hybes / reads.  
% 
% see also: 

defaults = cell(0,3);
defaults(end+1,:) = {'impute','boolean',false};
defaults(end+1,:) = {'removeEmpty','boolean',true};
defaults(end+1,:) = {'byHybe','boolean',true};
defaults(end+1,:) = {'dims','freeType',3};
defaults(end+1,:) = {'tableFormat',{'new','old','chrom'},'new'};
defaults(end+1,:) = {'chromCorrect','boolean',true};
defaults(end+1,:) = {'sortMethod',{'byPanel','byRead'},'byPanel'};
defaults(end+1,:) = {'dataBasedDriftCorrect','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

tableFormat = pars.tableFormat;

if isempty(spotTable)
    warning('spotTable is empty');
    distMap = [];
    polymer = [];
    return;
end

try

if pars.dataBasedDriftCorrect
    % remove fiduical based drift correction
    if quantile(abs(spotTable.xshift),.9) < 100
        % old version in pixels
       spotTable.x = spotTable.x + spotTable.xshift*154;
       spotTable.y = spotTable.y + spotTable.yshift*154;
       spotTable.z = spotTable.z + spotTable.zshift; 
    else
       spotTable.x = spotTable.x + spotTable.xshift;
       spotTable.y = spotTable.y + spotTable.yshift;
       spotTable.z = spotTable.z + spotTable.zshift;
    end
end
    
    
if strcmp(tableFormat,'chrom')
    uids = unique(spotTable.fs);
    numSpots = length(uids);
    numPanels = max(spotTable.panel);
    numReads = max(spotTable.readout);
    if strcmp(pars.sortMethod,'byPanel')
        distMap = nan(numPanels,numPanels,numSpots);
        polymer = nan(numPanels,4,numSpots);
    elseif strcmp(pars.sortMethod,'byRead')
        distMap = nan(numReads,numReads,numSpots);
        polymer = nan(numReads,4,numSpots);
    end
    
    
    spotProps = table(zeros(numSpots,1),zeros(numSpots,1),zeros(numSpots,1),...
        'VariableNames',{'locusX','locusY','fov'});
    for s=1:numSpots
        % isS1 = spotTable.s == s; % logical vector, all entries in table corresponding to spot s.      
        if strcmp(pars.sortMethod,'byPanel')
            isS1 = spotTable.fs == uids(s) ;% % the readout probe numbers of these spots
        elseif strcmp(pars.sortMethod,'byRead')
            isS1 = spotTable.fs == uids(s) & strcmp(spotTable.dataType,'H') ; % the readout probe numbers of these spot
        end      
        if pars.chromCorrect
            xyzh = [spotTable.x(isS1),spotTable.y(isS1),spotTable.z(isS1),spotTable.h(isS1)]  ...  % +   the drift-corrected (x,y,z) corrdinates
                -1*[spotTable.xcShift(isS1),spotTable.ycShift(isS1),spotTable.zcShift(isS1),0*spotTable.h(isS1)] ; % chromatic-correction       
        else
            xyzh = [spotTable.x(isS1),spotTable.y(isS1),spotTable.z(isS1),spotTable.h(isS1)];
        end  
        if strcmp(pars.sortMethod,'byRead')
            b = spotTable.readout(isS1);     
        elseif strcmp(pars.sortMethod,'byPanel')
            b = spotTable.panel(isS1);
        end
        if length(b)>1
            polymer(b,:,s) = xyzh;
            lx = spotTable.locusX(isS1);
            ly = spotTable.locusY(isS1);
            fov = spotTable.fov(isS1);
            spotProps(s,:) = {lx(1),ly(1),fov(1)};
            if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                sMap = squareform( pdist(xyzh(:,1:3)) );
            elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
                sMap = squareform( pdist(xyzh(:,1:2)) );
            elseif strcmp(pars.dims,'xz')
                sMap = squareform( pdist(xyzh(:,[1,3])) );
            elseif strcmp(pars.dims,'yz')
                sMap = squareform( pdist(xyzh(:,2:3)) );
            end
            distMap(b,b,s) = sMap;
        end
    end
    
elseif strcmp(tableFormat,'old')
    
    numSpots = max(spotTable.s);
    numReads = max(spotTable.hybe);
    distMap = nan(numReads,numReads,numSpots);
    polymer = nan(numReads,4,numSpots);
    spotProps = table(zeros(numReads,1),zeros(numReads,1),zeros(numReads,1),...
        'VariableNames',{'locusX','locusY','fov'});
    for s=1:numSpots
        isS1 = spotTable.s == s; % logical vector, all entries in table corresponding to spot s.
        xyzh = [spotTable.xc(isS1),spotTable.yc(isS1),spotTable.zc(isS1),spotTable.h(isS1)]; % the drift-corrected (x,y,z) corrdinates
        if pars.impute
           xyzh = fillmissing(xyzh,'linear'); 
        end
        b = spotTable.hybe(isS1); % the readout probe numbers of these spots
        polymer(b,:,s) = xyzh;
        lx = spotTable.locusX(isS1);
        ly = spotTable.locusY(isS1);
        fov = spotTable.fov(isS1);
        spotProps(s,:) = {lx(1),ly(1),fov(1)};      
        if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                sMap = squareform( pdist(xyzh(:,1:3)) );
        elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
            sMap = squareform( pdist(xyzh(:,1:2)) );
        elseif strcmp(pars.dims,'xz')
            sMap = squareform( pdist(xyzh(:,[1,3])) );
        elseif strcmp(pars.dims,'yz')
            sMap = squareform( pdist(xyzh(:,2:3)) );
        end
        
        distMap(b,b,s) = sMap;
    end

elseif strcmp(tableFormat,'new') && ~pars.byHybe  % new format  
    
    numSpots = max(spotTable.s);
    numReads = max(spotTable.read);
    distMap = nan(numReads,numReads,2*numSpots);
    polymer = nan(numReads,4,2*numSpots);
    spotProps = table(zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),...
        'VariableNames',{'locusX','locusY','fov','s'});
    
    for s=1:numSpots
        isS1 = spotTable.s == s & spotTable.idx == 1; % logical vector, all entries in table corresponding to spot s.
        isS2 = spotTable.s == s & spotTable.idx == 2; % logical vector, all entries in table corresponding to spot s.
        xyzh1 = [spotTable.x(isS1),spotTable.y(isS1),spotTable.z(isS1),spotTable.h(isS1)]; % the drift-corrected (x,y,z) corrdinates
        xyzh2 = [spotTable.x(isS2),spotTable.y(isS2),spotTable.z(isS2),spotTable.h(isS2)]; % the drift-corrected (x,y,z) corrdinates
        if pars.impute
           xyzh1 = fillmissing(xyzh1,'linear'); 
           xyzh2 = fillmissing(xyzh2,'linear'); 
        end
        s1 = 2*(s-1)+1;
        b = spotTable.read(isS1); % the readout probe numbers of these spots
        if length(b) > 1
            polymer(b,:,s1) = xyzh1;
  
            
            if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                sMap = squareform( pdist(xyzh1(:,1:3)) );
            elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
                sMap = squareform( pdist(xyzh1(:,1:2)) );
            elseif strcmp(pars.dims,'xz')
                sMap = squareform( pdist(xyzh1(:,[1,3])) );
            elseif strcmp(pars.dims,'yz')
                sMap = squareform( pdist(xyzh1(:,2:3)) );
            end
            
            
            distMap(b,b,s1) = sMap;
            lx = spotTable.locusX(isS1);
            ly = spotTable.locusY(isS1);
            fov = spotTable.fov(isS1);
            spotProps.locusX(s1) = lx(1);
            spotProps.locusY(s1) = ly(1);
            spotProps.fov(s1) =  fov(1);
            spotProps.s(s1) =  s;
        end
        s2 = 2*s;
        b = spotTable.read(isS2); % the readout probe numbers of these spots
        if length(b) > 1
            polymer(b,:,s2) = xyzh2;
            lx = spotTable.locusX(isS2);
            ly = spotTable.locusY(isS2);
            fov = spotTable.fov(isS2);
            spotProps.locusX(s2) = lx(1);
            spotProps.locusY(s2) = ly(1);
            spotProps.fov(s2) =  fov(1);
            spotProps.s(s2) =  s;
            
            if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                sMap = squareform( pdist(xyzh2(:,1:3)) );
            elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
                sMap = squareform( pdist(xyzh2(:,1:2)) );
            elseif strcmp(pars.dims,'xz')
                sMap = squareform( pdist(xyzh2(:,[1,3])) );
            elseif strcmp(pars.dims,'yz')
                sMap = squareform( pdist(xyzh2(:,2:3)) );
            end
            
            
            if ~isempty(sMap)
                distMap(b,b,s2) = sMap;
            end
        end
    end
    if pars.removeEmpty
        blank = squeeze(nansum(polymer(:,1,:),1)) == 0;
        distMap(:,:,blank) = [];
        polymer(:,:,blank) = [];
        spotProps(blank,:) = [];
    end
    
elseif pars.byHybe    
    
    numSpots = max(spotTable.s);
    numHybs = max(spotTable.hybe);
    distMap = nan(numHybs,numHybs,2*numSpots);
    polymer = nan(numHybs,4,2*numSpots);
    spotProps = table(zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),...
        'VariableNames',{'locusX','locusY','fov','s'});
    for s=1:numSpots
        isS1 = spotTable.s == s & spotTable.idx == 1; % logical vector, all entries in table corresponding to spot s.
        isS2 = spotTable.s == s & spotTable.idx == 2; % logical vector, all entries in table corresponding to spot s.
        xyzh1 = [spotTable.x(isS1),spotTable.y(isS1),spotTable.z(isS1),spotTable.h(isS1)]; % the drift-corrected (x,y,z) corrdinates
        xyzh2 = [spotTable.x(isS2),spotTable.y(isS2),spotTable.z(isS2),spotTable.h(isS2)]; % the drift-corrected (x,y,z) corrdinates
        if pars.impute
            error('imputation code is wrong. spots need to be sorted into their bins first!')
           xyzh1 = fillmissing(xyzh1,'linear'); 
           xyzh2 = fillmissing(xyzh2,'linear'); 
        end
        s1 = 2*(s-1)+1;
        b = spotTable.hybe(isS1); % the readout probe numbers of these spots
        if length(b) > 1
            polymer(b,:,s1) = xyzh1;
            lx = spotTable.locusX(isS1);
            ly = spotTable.locusY(isS1);
            fov = spotTable.fov(isS1);
            spotProps.locusX(s1) = lx(1);
            spotProps.locusY(s1) = ly(1);
            spotProps.fov(s1) =  fov(1);
            spotProps.s(s1) =  s;
           %  sMap = squareform( pdist(xyzh1(:,1:3)) );
          
            if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                    sMap = squareform( pdist(xyzh1(:,1:3)) );
            elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
                sMap = squareform( pdist(xyzh1(:,1:2)) );
            elseif strcmp(pars.dims,'xz')
                sMap = squareform( pdist(xyzh1(:,[1,3])) );
            elseif strcmp(pars.dims,'yz')
                sMap = squareform( pdist(xyzh1(:,2:3)) );
            end
            
            
            distMap(b,b,s1) = sMap;
        end
        s2 = 2*s;
        b = spotTable.hybe(isS2); % the readout probe numbers of these spots
        if length(b) > 1
            polymer(b,:,s2) = xyzh2;
            lx = spotTable.locusX(isS2);
            ly = spotTable.locusY(isS2);
            fov = spotTable.fov(isS2);
            spotProps.locusX(s2) = lx(1);
            spotProps.locusY(s2) = ly(1);
            spotProps.fov(s2) =  fov(1);
            spotProps.s(s2) =  s;
            if all(pars.dims == 3) || strcmp(pars.dims,'xyz')
                sMap = squareform( pdist(xyzh2(:,1:3)) );
            elseif all(pars.dims == 2) || strcmp(pars.dims,'xy')
                sMap = squareform( pdist(xyzh2(:,1:2)) );
            elseif strcmp(pars.dims,'xz')
                sMap = squareform( pdist(xyzh2(:,[1,3])) );
            elseif strcmp(pars.dims,'yz')
                sMap = squareform( pdist(xyzh2(:,2:3)) );
            end
            
            if ~isempty(sMap)
                distMap(b,b,s2) = sMap;
            end
        end
    end
    spotProps.uniqueID = spotProps.s + spotProps.fov*1E5;
    if pars.removeEmpty
        blank = squeeze(nansum(polymer(:,1,:),1)) == 0;
        distMap(:,:,blank) = [];
        polymer(:,:,blank) = [];
        spotProps(blank,:) = [];
    end
    
end

if pars.dataBasedDriftCorrect
    % compute typical offset
    [numHybes,~,nPolys] = size(polymer);
    dPoly = nan(numHybes,3,nPolys);
    for h=1:numHybes
       dPoly(h,:,:) = polymer(h,1:3,:)-polymer(1,1:3,:); 
    end
    meanOff = nanmedian(dPoly(:,1:3,:),3);
    % new polymer
    newPolymer = polymer;
    newPolymer(:,1:3,:)=newPolymer(:,1:3,:) - repmat(meanOff,[1,1,nPolys]);
    newMap = distMap;
    for n=1:nPolys
        if pars.dims == 3
            newMap(:,:,n) = squareform( pdist( newPolymer(:,1:3,n)  ));
        else
            newMap(:,:,n) = squareform( pdist( newPolymer(:,1:2,n)  ));
        end
    end
    % export
    polymer = newPolymer;
    distMap = newMap;
end



catch er
    error(er.getReport);
end
    