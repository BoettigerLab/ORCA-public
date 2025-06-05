function [Xs,Ys,Zs,Vs] = PlotLEF(xyz1,xyz2,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'maxSep','nonnegative',2.5};
defaults(end+1,:) = {'LEFweight','nonnegative',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);


%%
% N = 5;
% s = randi(10,N,3);
% xyz1 = s++rand(N,3);
% xyz2 = s+2+rand(N,3);

N = size(xyz1,1);
X = cell(N,1);
Y = cell(N,1);
Z = cell(N,1);
V = cell(N,1);

for r=1:N
    v = xyz2(r,:)-xyz1(r,:);
    d = sqrt( sum((xyz1(r,:)-xyz2(r,:)).^2) );
    if d<pars.maxSep 
        ring = [xyz1(r,:) - .25*v; 
                (xyz1(r,:) + .5*v.*[1,1,1/2])  ;  % (xyz1(r,:) + .5*v.*[1,1,1/2])  ;            ;
                xyz2(r,:) + .25*v;
                (xyz2(r,:) - .5*v.*[1,1,-1/2]);  % (xyz2(r,:) - .5*v.*[1,1,1/2])  ;            
                xyz1(r,:) - .25*v];
        
        [X{r},Y{r},Z{r},V{r}] = PlotTube(ring,'r',pars.LEFweight*.2*d,'interpPts',10,'method','spline','colormap',[1 .5 0],'plot',false);
        if r<N
            X{r} = [X{r}; nan(1,size(X{r},2))];
            Y{r} = [Y{r}; nan(1,size(X{r},2))];
            Z{r} = [Z{r}; nan(1,size(X{r},2))];
            V{r} = [V{r}; nan(1,size(X{r},2))];
        end
    end
end
Xs = cat(1,X{:});
Ys = cat(1,Y{:});
Zs = cat(1,Z{:});
Vs = cat(1,V{:});

if nargout == 0
    surf(Xs,Ys,Zs,Vs);
end


