function mRNAwarp = ComputeWarpByMRNA(cents,datPos,xf,yf,varargin)
%

% Mcents

% default Parameters
showPlots = true; 
maxDtoCentroid = .5; % max distance you can be from the centroid of the mRNA
goodMRNA = [];

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 4
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'maxDtoCentroid'
                maxDtoCentroid = CheckParameter(parameterValue,'positive','maxDtoCentroid');
            case 'goodMRNA'
                goodMRNA = CheckParameter(parameterValue,'cell','goodMRNA');
            case 'showPlots'
                showPlots = CheckParameter(parameterValue,'boolean','showPlots');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%% Main Function 
if showPlots
    figure(6); clf;
    plot(0,0,'ko'); hold on;
    plot(0,0,'k.');
    plot(0,0,'k+');
    plot(cents(:,1),cents(:,2),'k.'); hold on;
    figure(4); clf;
end

tform_start = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);
method = 'nonreflective similarity';

D = length(datPos); 
cMap = jet(D); 
mRNAwarp = cell(D,1); 
for d=1:D
    if isempty(goodMRNA)
        refPos = cents; % This works alright
    else
        refPos = cents(goodMRNA{d},:);
    end
   %   refPos = cents(goodMRNA,:);
   
    Ref.x = refPos(:,1);
    Ref.y = refPos(:,2); 
    Dat.x = datPos{d}(:,1);
    Dat.y = datPos{d}(:,2);
    
    [matched, ~] = corr_mols(Ref, Dat, tform_start, maxDtoCentroid);
    RefP = [Ref.x(matched.set1_inds),Ref.y(matched.set1_inds)];
    DatP = [Dat.x(matched.set2_inds),Dat.y(matched.set2_inds)];   

%         figure(5); clf;
%         plot(cents(:,1),cents(:,2),'k.'); hold on;
%         plot(xf{d},yf{d},'r.','MarkerSize',1); hold on;
%         plot(Ref.x,Ref.y,'k+'); hold on;
%         plot(Dat.x,Dat.y,'r.'); hold on;
%         plot(RefP(:,1),RefP(:,2),'ko','MarkerSize',10);
%         plot(DatP(:,1),DatP(:,2),'ro','MarkerSize',10);
    
     mRNAwarp{d} = cp2tform(RefP,DatP,method); % compute warp
    % mRNAwarp{d} = cp2tform(RefP,DatP,'polynomial',3); % compute warp
    
    if showPlots
        [xw,yw] = tforminv(mRNAwarp{d},DatP(:,1),DatP(:,2));
        figure(6); 
        plot(Ref.x,Ref.y,'k+'); hold on;
        plot(Dat.x,Dat.y,'color',cMap(d,:),'Marker','o','LineStyle','none','MarkerSize',10); hold on;
        plot(RefP(:,1),RefP(:,2),'color','k','Marker','+','LineStyle','none','MarkerSize',5); hold on;
        plot(DatP(:,1),DatP(:,2),'color',cMap(d,:),'Marker','o','LineStyle','none','MarkerSize',10); hold on;
        plot(xw,yw,'color',cMap(d,:),'Marker','.','LineStyle','none','MarkerSize',10); hold on;
        plot(xf{d},yf{d},'color',cMap(d,:),'Marker','x','LineStyle','none','MarkerSize',4); hold on;
        
        N = length(xw);
        links = zeros(3*N,2);
        links(1:3:end,:) = [xw,yw];
        links(2:3:end,:) = RefP;
        links(3:3:end,:) = NaN*ones(N,2);
        plot(links(:,1),links(:,2),'k');
        
        N = length(xw);
        links1 = zeros(3*N,2);
        links1(1:3:end,:) = DatP;
        links1(2:3:end,:) = RefP;
        links1(3:3:end,:) = NaN*ones(N,2);
        plot(links1(:,1),links1(:,2),'r');
        
        
        
        figure(5); clf; quiver(DatP(:,1),DatP(:,2),xw-DatP(:,1),yw-DatP(:,2));
        
        orig_err = sqrt( (RefP(:,1)-DatP(:,1)).^2 + (RefP(:,2)-DatP(:,2)).^2 );
        warp_err = sqrt( (xw-DatP(:,1)).^2 + (yw-DatP(:,2)).^2 );     
        
        figure(4); subplot(round(D/2),2,d);
        hist(orig_err,linspace(0,maxDtoCentroid,100)); hold on;
        hist(warp_err,linspace(0,maxDtoCentroid,100));
        hs = findobj('Type','patch');
        set(hs(1),'FaceColor','r','EdgeColor','r'); 
    end   
end

if showPlots
    figure(6);     
    legend({'original spot','warped spot','mRNA centroid'});
end

% newWarp = cell(D,1); 
% for d=1:D
%    newWarp{d} = maketform('composite',tform{d},mRNAwarp{d}) ;
% end
%      
