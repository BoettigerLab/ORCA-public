function newMap = InterpMapNans(map,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'badHybes','freeType',[]};
defaults(end+1,:) = {'badPixels','freeType',[]};
defaults(end+1,:) = {'method',{'nextNonempty','moving','interp','loess','rloess'},'moving'}; % nextNonempty
defaults(end+1,:) = {'window','nonnegative',4}; % nextNonempty
defaults(end+1,:) = {'smoothMethod',{'moving','lowess','loess','sgolay','rlowess','rloess'},'moving'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if ~isempty(pars.badHybes)
    map(pars.badHybes,:) = NaN;
    map(:,pars.badHybes) = NaN;
end
if ~isempty(pars.badPixels)
    for p=1:size(pars.badPixels,1)
        map(pars.badPixels(p,1),pars.badPixels(p,2)) = NaN;
        map(pars.badPixels(p,2),pars.badPixels(p,1)) = NaN;
    end
end
% figure(10); clf; imagesc(map);
nReads = size(map,1);

% find a nonzero row, use it to get all the nan columns
n = 1; 
stop = false;
while n<nReads && ~stop
    if nansum(map(n,:)) <= 1
        n=n+1;
    else
        stop = true;
    end
end
if ~isempty(pars.badPixels)
    badCols = [find(isnan(map(n,:))),pars.badPixels(:,1)',pars.badPixels(:,2)'];
else
    badCols = find(isnan(map(n,:)));
end
% don't interpolate if there isn't data on either side
nRows = nReads;
badEnd = ismember(badCols,nRows);
while any(badEnd)
    badEnd = ismember(badCols,nRows);
    badCols(badEnd) = [];
    nRows = nRows -1;
end

newMap = map;
tempMap = map;
% figure(10); clf; imagesc(newMap);

if strcmp(pars.method,'moving') || strcmp(pars.method,'loess') || strcmp(pars.method,'rloess') 
    
    for r=1:nRows
        yy = map(r,1:nRows);
        if any(~isnan(yy))
             yvalues= smooth(1:nRows,yy,pars.window,pars.method);
            tempMap(r,1:nRows) = yvalues;
        end
    end
    % figure(10); clf; imagesc(tempMap);
    for r=1:nRows
        yy = tempMap(1:nRows,r);
        if any(~isnan(yy))
            yvalues= smooth(1:nRows,yy,pars.window,pars.method);
            tempMap(1:nRows,r) = yvalues;
        end
    end
%     figure(10); clf; imagesc(tempMap);
    
    for r=badCols
        newMap(1:nRows,r) = tempMap(1:nRows,r);
        newMap(r,1:nRows) = tempMap(r,1:nRows);
    end
%     figure(10); clf;
%     s1 = subplot(1,2,1); imagesc(map);
%     s2 = subplot(1,2,2); imagesc(newMap);

elseif strcmp(pars.method,'interp')
    
    for r=1:nRows
        yy = map(r,1:nRows);
        x = 1:nRows;
        xx = x; 
        xx(badCols) = [];
        yy(badCols) = [];
        yvalues = interp1(xx,yy,x);
        tempMap(r,1:nRows) = yvalues;
    end
    % figure(10); clf; imagesc(tempMap);
    for r=1:nRows
        yy = tempMap(1:nRows,r);
        x = 1:nRows;
        xx = x; 
        xx(badCols) = [];
        yy(badCols) = [];
        yvalues = interp1(xx,yy,x);
        tempMap(1:nRows,r) = yvalues;
    end
%     figure(10); clf; imagesc(tempMap);
    
    for r=badCols
        newMap(1:nRows,r) = tempMap(1:nRows,r);
        newMap(r,1:nRows) = tempMap(r,1:nRows);
    end
%     figure(10); clf;
%     s1 = subplot(1,2,1); imagesc(map);
%     s2 = subplot(1,2,2); imagesc(newMap);
    
elseif strcmp(pars.method,'nextNonempty')
    for b=badCols
        % grab nonzero columns to the left (if there is a left)
        n = b-1;
        stop = false;
        leftVals = map(:,b);
        while n > 0 && ~stop   
           if nansum(map(:,n)) == 0
               n=n-1;
           else
               stop = true;
               leftVals = map(:,n);
           end       
        end

        n = b+1;
        stop = false;
        rightVals = map(:,b);
        while n < nReads && ~stop   
           if nansum(map(:,n)) == 0
               n=n+1;
           else
               stop = true;
               rightVals = map(:,n);
           end
        end

        newVals =  nanmean([leftVals,rightVals],2);      
        newMap(:,b) = newVals;
        % newMap(b,:) = newVals;
    end


    % loop over rows
    badRows = badCols;
    for b=badRows
        % grab nonzero columns to the left (if there is a left)
        n = b-1;
        stop = false;
        leftVals = newMap(b,:);
        while n > 0 && ~stop   
           if nansum(newMap(n,:)) == 0
               n=n-1;
           else
               stop = true;
               leftVals = newMap(n,:);
           end       
        end

        n = b+1;
        stop = false;
        rightVals = newMap(b,:);
        while n < nReads && ~stop   
           if nansum(newMap(n,:)) == 0
               n=n+1;
           else
               stop = true;
               rightVals = newMap(n,:);
           end
        end

        newVals =  nanmean([leftVals;rightVals],1);      
        newMap(b,:) = newVals;

    % figure(1); clf; imagesc(newMap); pause(.1);
    end
% figure(2); clf; imagesc(newMap);
end
