function [distMap,polymer,spotProps] = TableToDistMap(spotTable,varargin)
% distMap = TableToDistMap(bxcTable)
% sorts by hybes / reads.  

defaults = cell(0,3);
defaults(end+1,:) = {'impute','boolean',false};
defaults(end+1,:) = {'removeEmpty','boolean',true};
defaults(end+1,:) = {'byHybe','boolean',true};
defaults(end+1,:) = {'dims','integer',3};

pars = ParseVariableArguments(varargin,defaults,mfilename); 

try
    spotTable.xc(1);
    tableFormat = 'old';
catch
    tableFormat = 'new';
end

if isempty(spotTable)
    warning('spotTable is empty');
    distMap = [];
    polymer = [];
    return;
end

try

if strcmp(tableFormat,'old')
    
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
%         spotProps.locusX(s) = lx(1);
%         spotProps.locusY(s) = ly(1);
%         spotProps.fov(s) =  fov(1);
        if pars.dims  == 3
            sMap = squareform( pdist(xyzh(:,1:3)) );
        elseif pars.dims == 2
            sMap = squareform( pdist(xyzh(:,1:2)) );
        end
        distMap(b,b,s) = sMap;
    end

elseif strcmp(tableFormat,'new') && ~pars.byHybe  % new format  
    
    numSpots = max(spotTable.s);
    numReads = max(spotTable.read);
    distMap = nan(numReads,numReads,2*numSpots);
    polymer = nan(numReads,4,2*numSpots);
    spotProps = table(zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),...
        'VariableNames',{'locusX','locusY','fov'});
    
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
  
            if pars.dims  == 3
                sMap = squareform( pdist(xyzh1(:,1:3)) );
            elseif pars.dims == 2
                sMap = squareform( pdist(xyzh1(:,1:2)) );
            end
            
            distMap(b,b,s1) = sMap;
            lx = spotTable.locusX(isS1);
            ly = spotTable.locusY(isS1);
            fov = spotTable.fov(isS1);
            spotProps.locusX(s1) = lx(1);
            spotProps.locusY(s1) = ly(1);
            spotProps.fov(s1) =  fov(1);
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
            
             if pars.dims  == 3
                sMap = squareform( pdist(xyzh2(:,1:3)) );
            elseif pars.dims == 2
                sMap = squareform( pdist(xyzh2(:,1:2)) );
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
    spotProps = table(zeros(2*numSpots,1),zeros(2*numSpots,1),zeros(2*numSpots,1),...
        'VariableNames',{'locusX','locusY','fov'});
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
        b = spotTable.hybe(isS1); % the readout probe numbers of these spots
        if length(b) > 1
            polymer(b,:,s1) = xyzh1;
            lx = spotTable.locusX(isS1);
            ly = spotTable.locusY(isS1);
            fov = spotTable.fov(isS1);
            spotProps.locusX(s1) = lx(1);
            spotProps.locusY(s1) = ly(1);
            spotProps.fov(s1) =  fov(1);
           %  sMap = squareform( pdist(xyzh1(:,1:3)) );
            
            if pars.dims  == 3
                sMap = squareform( pdist(xyzh1(:,1:3)) );
            elseif pars.dims == 2
                sMap = squareform( pdist(xyzh1(:,1:2)) );
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
            
           %  sMap = squareform( pdist(xyzh2(:,1:3)) );
            if pars.dims  == 3
                sMap = squareform( pdist(xyzh2(:,1:3)) );
            elseif pars.dims == 2
                sMap = squareform( pdist(xyzh2(:,1:2)) );
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
    
end

catch er
    error(er.getReport);
end
    