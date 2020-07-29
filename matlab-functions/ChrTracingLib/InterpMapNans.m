function newMap = InterpMapNans(map,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'badHybes','freeType',[]};
defaults(end+1,:) = {'method',{'nextNonempty','interp'},'nextNonempty'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if ~isempty(pars.badHybes)
    map(pars.badHybes,:) = NaN;
    map(:,pars.badHybes) = NaN;
end

nReads = size(map,1);

% find a nonzero row, use it to get all the nan columns
n = 1; 
stop = false;
while n<nReads && ~stop
    if nansum(map(n,:)) == 0
        n=n+1;
    else
        stop = true;
    end
end
badCols = find(isnan(map(n,:)));

newMap = map;

if strcmp(pars.method,'interp')
    % not working yet: use a longer-range averaging filter
    noData = isnan(map);
    nRows = size(map,1);
    for r=1:nRows
        if ~any(r==badCols)
            yy = map(r,:);
            yvalues= smooth(yy,0.05,'rloess');
            addData = yvalues;
            addData(~noData(r,:)) = map(r,~noData(r,:))';
            % xx = xvalues;
            % yy(noData(r,:)) = [];
            % xx(noData(r,:)) = [];
            % yvalues= smooth(yy,0.05,'rloess');
            % yvalues = interp1(xx,yvalues,xvalues,'spline');
            newMap(r,:) = addData;
            % figure(3); clf; plot(xx,log10(yy),'k.');
            % hold on; plot(xvalues,log10(yvalues),'ro');
            % pause;
        end
    end
    figure(3); clf; subplot(1,2,1); imagesc(log10(map));
    newMap(newMap<0.001) = 0;
    subplot(1,2,2); imagesc(log10(newMap));
    
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
