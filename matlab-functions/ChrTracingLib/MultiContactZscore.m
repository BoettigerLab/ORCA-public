function zval = MultiContactZscore(polstack,mapstack,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'threshold','positive',200};
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);
%%

[n,~,nTraces] = size(polstack);

disBA = zeros(n,n,n);
disBC = zeros(n,n,n);
for a=1:n
    for b=1:n
        for c=1:n
            disBA(a,b,c) = abs(a-b);
            disBC(a,b,c) = abs(c-b);
        end
    end
end

%%
% this will be slow
theta = pars.threshold; % cutoff in nanometers
C = size(polstack,3);
contactCube = zeros(n,n,n);
for a=1:n
    if pars.verbose
    disp([num2str(a/n*100),'% complete']);
    end
    for b=1:n
        for c=1:n % a=4, b=4, c =7 
            contactCube(a,b,c) = sum(mapstack(a,b,:) < theta & mapstack(c,b,:) < theta) ./ (sum(mapstack(a,b,:) < inf & mapstack(c,b,:) < inf)  ); 
        end
    end
end

% a = 3;
% raw3 = squeeze(contactCube(a,:,:));
% figure(1); clf; imagesc(raw3); colorbar; caxis([0,.1]);
%% Normalize
if pars.verbose
disp('normalizing data...')
end
% loop over different distance intervals to select all distance matched comparison sites 
zval = zeros(n,n,n);
for a=1:n
    if pars.verbose
        disp(['normalizing data... ',num2str(a/n*100),'% complete']);
    end
    for x=1:n-1
        for y=1:n-1 % a=4; x = 4; y=7;
             match = disBA==abs(a-x) & disBC==abs(a-y);
             compPt = contactCube(match);
             zval(a,x,y) = (contactCube(a,x,y) - nanmean(compPt))/nanstd(compPt);
        end
    end
end

%%
% a = 15;
% figure(1); clf; 
% im = squeeze(contactCube(a,:,:)); %  im(im>.05) = NaN; 
% subplot(1,2,1); imagesc(im); colorbar; caxis([0,.08]);
% subplot(1,2,2); imagesc(squeeze(zval(a,:,:))); colorbar;   caxis([-4,4]);
% GetColorMap('BlueWhiteRedSat')