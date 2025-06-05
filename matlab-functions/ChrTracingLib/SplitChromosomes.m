function multiFits2 = SplitChromosomes(multiFits1,varargin)
% Sort two polymer traces by swapping monomers so that both have the
% shortest path and maximal separation possible. 
% starts at the first non-missing data point and adds to minimize the sum
% of the distance of both growing paths. 
% 
%% Inputs
% multifits = N-monomers x 3-dims x 2-polymers matrix
% 
%% Outputs
% multifits = N-monomers x 3-dims x 2-polymers matrix

defaults = cell(0,3);
defaults(end+1,:) = {'maxStep','positive',inf}; % steps beyond this distance will be treated as NaNs and ignored to avoid confusing the chain
defaults(end+1,:) = {'showPlots','boolean',false}; % show comparison of initial and final polymers
defaults(end+1,:) = {'showExtraPlots','boolean',false}; % show progress as spot pairs are assigned to separate traces / polymers 
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% function [polymer1Out,polymer2Out] = SplitChromosomes(polymer1,polymer2,varargin)
% multiFits1 = cat(3,polymer1,polymer2);


% build two synthetically separate polys
% multiFits1 = multiFits;
% multiFits1(:,:,1) = multiFits1(:,:,1)+8000;
% b = rand(nB,1)>.5;
% multiFitsTemp = multiFits1;
% multiFits1(b,:,1) = multiFitsTemp(b,:,2);
% multiFits1(b,:,2) = multiFitsTemp(b,:,1);

nB = size(multiFits1,1); 
if pars.verbose
    % compute original path difference
    m12 = squareform(pdist(cat(1,multiFits1(:,:,1),multiFits1(:,:,2))));
    [nanmean(m12(1:nB,1:nB),'all'),nanmean(m12(nB+1:2*nB,nB+1:2*nB),'all'),nanmean(m12(1:nB,nB+1:2*nB),'all')]
    sep1 =nanmean(m12(1:nB,nB+1:2*nB),'all');
end

% figure(6); clf; imagesc(m12); colorbar; clim([0,1000])


multiFits2=multiFits1;
a=1; % current index
b=2; % index of next pair to assign into paths
if pars.showExtraPlots
    figure(5); clf; 
    plot(multiFits2(a,1,1),multiFits2(a,2,1),'r.'); hold on;
    plot(multiFits2(a,1,2),multiFits2(a,2,2),'b.'); hold on;
    text(multiFits2(a,1,1),multiFits2(a,2,1),num2str(1),'color','r');
    text(multiFits2(a,1,2),multiFits2(a,2,2),num2str(1),'color','b');
    plot(multiFits2(b,1,1),multiFits2(b,2,1),'r.'); hold on;
    plot(multiFits2(b,1,2),multiFits2(b,2,2),'b.'); hold on;
    text(multiFits2(b,1,1),multiFits2(b,2,1),num2str(2),'color','r');
    text(multiFits2(b,1,2),multiFits2(b,2,2),num2str(2),'color','b');
end
% init xyz0 should probably be the centroid. 
while b<nB % b=1 
    xyz0 = squeeze(multiFits2(a,:,1:2))'; 
    xyz1 = squeeze(multiFits2(b,:,1:2))';
     if any(isnan(xyz0(:))) % should only kick in if we don't have two spots for the initial point  
         a=a+1;
         b=b+1;
     else
        if any(isnan(xyz1(:))) % if we don't have data to match keep going 
            b=b+1;
        else
            [m,cost] = MatchPoints(xyz0,xyz1);
            if sum(cost) > pars.maxStep % 2*3e3 % max step
                % we want to ignore this point, so we go to the next and
                % leave the source (a) where it was
                b=b+1;
            else
                % we make the swap
                if pars.showExtraPlots
                    figure(5); 
                    plot(multiFits2(b,1,1),multiFits2(b,2,1),'ro'); hold on;
                    plot(multiFits2(b,1,2),multiFits2(b,2,2),'bo'); hold on;
                end
                multiFits2(b,:,m(:,2))=multiFits2(b,:,m(:,1));
                if m(1,1)==2 && pars.verbose
                    disp(['swapped b=',num2str(b)]);
                end
                if pars.showExtraPlots
                    plot(multiFits2(b,1,1),multiFits2(b,2,1),'r.'); hold on;
                    plot(multiFits2(b,1,2),multiFits2(b,2,2),'b.'); hold on;
                    text(multiFits2(b,1,1),multiFits2(b,2,1),num2str(b),'color','r');
                    text(multiFits2(b,1,2),multiFits2(b,2,2),num2str(b),'color','b');
                    pause(.1);
                end
                a=b; % we update this point to be the new source
                b=b+1; % we increase the next point
            end
        end
     end
end

if pars.verbose
    % compute improvement
    m12 = squareform(pdist(cat(1,multiFits2(:,:,1),multiFits2(:,:,2))));
    [nanmean(m12(1:nB,1:nB),'all'),nanmean(m12(nB+1:2*nB,nB+1:2*nB),'all'),nanmean(m12(1:nB,nB+1:2*nB),'all')]
    sep2 =nanmean(m12(1:nB,nB+1:2*nB),'all');
    % display improvement
    disp(['more separated by ',num2str(100*(1-sep1/sep2),5),'%']);
end

% polymer1Out = multiFits2(:,:,1);
% polymer2Out = multiFits2(:,:,2); 

if pars.showPlots
    figure(7); clf; 
    subplot(1,2,1);
    plot3(multiFits1(:,1,1),multiFits1(:,2,1),multiFits1(:,3,1),'.-'); hold on;
    plot3(multiFits1(:,1,2),multiFits1(:,2,2),multiFits1(:,3,2),'.-'); hold on;
    view([90,0])
    subplot(1,2,2);
    plot3(multiFits2(:,1,1),multiFits2(:,2,1),multiFits2(:,3,1),'.-'); hold on;
    plot3(multiFits2(:,1,2),multiFits2(:,2,2),multiFits2(:,3,2),'.-'); hold on;
    view([90,0])
end