function tform = ComputeWarpByBeads(FedPos,xshift,yshift,varargin)

showPlots = true; 

D = length(xshift); % number of stains 
cMap = jet(D);

tform_start = maketform('affine',[1 0 0; 0 1 0; 0 0 1]);
match_radius =20;
method = 'nonreflective similarity';
tform = cell(D,1); 
Ref.x = xshift(1) + FedPos{1}(:,1);
Ref.y = yshift(1) + FedPos{1}(:,2);

% Compute distances to all neighbors.
% Keep beads which have the same inter-neigbor distance relationship in
% each frame.  Uses these beads for the warp map. 
% Nbeads = length(Ref.x)

figure(3); clf; figure(4); clf;
for d=1:D % d=2
    Dat.x = xshift(d) + FedPos{d}(:,1);
    Dat.y = yshift(d) + FedPos{d}(:,2);
    [matched, ~] = corr_mols(Ref, Dat,tform_start, match_radius);
    RefP = [Ref.x(matched.set1_inds),Ref.y(matched.set1_inds)];
    DatP = [Dat.x(matched.set2_inds),Dat.y(matched.set2_inds)];
    
    if showPlots
        figure(3); 
        tform{d} = cp2tform(RefP,DatP,method); % compute warp
        plot(DatP(:,1),DatP(:,2),'color',cMap(d,:),'Marker','o','LineStyle','none','MarkerSize',20); hold on;
        [xw,yw] = tforminv(tform{d},DatP(:,1),DatP(:,2));
        plot(xw,yw,'color',cMap(d,:),'Marker','.','LineStyle','none','MarkerSize',50/d); hold on;

        orig_err = sqrt( (RefP(:,1)-DatP(:,1)).^2 + (RefP(:,2)-DatP(:,2)).^2 );
        warp_err = sqrt( (xw-DatP(:,1)).^2 + (yw-DatP(:,2)).^2 );
        figure(4); subplot(round(D/2),2,d);
        hist(orig_err,linspace(0,match_radius,10)); hold on;
        hist(warp_err,linspace(0,match_radius,10));
        hs = findobj('Type','patch');
        set(hs(1),'FaceColor','r');
    end
end