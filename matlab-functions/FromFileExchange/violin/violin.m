%--------------------------------------------------------------------------
%violin.m - Simple violin plot using matlab default kernel density estimation
%v2: extended for accepting also cells of different length
%--------------------------------------------------------------------------
%This function creates violin plots based on kernel density estimation
%using ksdensity with default settings. Please be careful when comparing pdfs 
%estimated with different bandwidth!
%--------------------------------------------------------------------------
% Modified from:
%Please cite this function as:
%Hoffmann H, 2013: violin.m - Simple violin plot using matlab default kernel 
%density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
%hhoffmann@uni-bonn.de
%--------------------------------------------------------------------------
%Input:
%X:     either:
% n x m matrix. A 'violin' is plotted for each column m, OR
% 1 x m Cellarry with elements being numerical colums of nx1 length.
%xL:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
%
%varargin:
%fc=[1 0.5 0]%FaceColor: Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
%lc='k'      %LineColor: Specify abbrev. (e.g. 'k' for black)
%alp=0.5     %Alpha value (transparency)   
%mc='k'      %Color of the bars indicating the mean
%medc='r'    %Color of the bars indicating the median
%--------------------------------------------------------------------------
%{
%Example1 (default):
disp('this example uses the statistical toolbox')
X=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
[h,L,MX,MED]=violin(X,[]); 
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%
%Example2 (specify):
disp('this example uses the statistical toolbox')
X=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
[h,L,MX,MED]=violin(X,{'a','b','c','d'},[1 1 0;0 1 0;.3 .3 .3;0 0.3 1],'w',0.3,'k','r--')
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%close all, [h,L,MX,MED]=violin(X,{'a','b','c','d'},[0.1 0.3 0.3],'k',1,'w','k--')
%}
%--------------------------------------------------------------------------
function[h,L,MX,MED]=violin(X,varargin)
 %defaults:

 
 defaults = cell(0,3);
defaults(end+1,:) = {'xOffset','integer',0};
defaults(end+1,:) = {'lineColor','freeType','k'};
defaults(end+1,:) = {'faceColor','freeType',[.5 .5 .5]};
defaults(end+1,:) = {'alpha','fraction',1};
defaults(end+1,:) = {'meanColor','freeType','k'};
defaults(end+1,:) = {'medColor','freeType','r'};
defaults(end+1,:) = {'plotMean','boolean',true};
defaults(end+1,:) = {'medianStyle','cell',{}};
defaults(end+1,:) = {'variableNames','cell',{''}};
defaults(end+1,:) = {'nonnegative','boolean',false};
defaults(end+1,:) = {'method',{'linear','nearest','next','pchip','cubic','makima','spline'},'linear'};
defaults(end+1,:) = {'bandwidth','positive',.05}; % on the scale of the data
pars = ParseVariableArguments(varargin,defaults,mfilename);
 

nData = size(X,2);
 lc= pars.lineColor; % 'k';
 fc= pars.faceColor; % [1 0.5 0];
 alp= pars.alpha; % 0.5;
 mc= pars.meanColor;%  'k';
 medc= pars.medColor;% 'r';
 xL = pars.variableNames;
 if isempty(xL)
     xL = cellstr(num2str( (1:nData)' ));
 end
 
 %convert everything to cells:
 if iscell(X)==0
     i=1;
     for i=1:nData
        X2{i}=X(:,i);
     end
      X=X2;
 end
 
 for n=1:nData
    temp = X{n};
    temp(isnan(temp)) = [];
    X{n} = temp;
 end
 
 %% Main Function
 
 %-------------------------------------------------------------------------
 if size(fc,1)==1
     fc=repmat(fc,size(X,2),1);
 end
 %-------------------------------------------------------------------------
 i=1;
 for i=1:size(X,2)
    [f, u]=ksdensity(X{i},'bandwidth',pars.bandwidth*quantile(X{i},.99));
    f=f/max(f)*0.3; %   normalize width
    F(:,i)=f;
    U(:,i)= u;
    MED(:,i)=nanmedian(X{i});
    MX(:,i)=nanmean(X{i});
 end

 %-------------------------------------------------------------------------
if pars.nonnegative
    U(U<0) = min(U(U>0));
end


 i=1;
 medError = false;
 for i=i:size(X,2)
    fill(pars.xOffset + [F(:,i)+0.5*i*2;flipud(2*i*0.5-F(:,i))],...
        [U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',alp,'EdgeColor',lc)
    hold on
    if pars.plotMean
        p(1)=plot(pars.xOffset + [interp1(U(:,i),F(:,i)+0.5*i*2,MX(:,i),pars.method),...
            interp1(flipud(U(:,i)),flipud(2*i*0.5-F(:,i)),MX(:,i),pars.method) ],...
            [MX(:,i) MX(:,i)],mc,'LineWidth',2);

        try
        p(2)=plot(pars.xOffset + [interp1(U(:,i),F(:,i)+0.5*i*2,MED(:,i)),...
            interp1(flipud(U(:,i)),flipud(2*i*0.5-F(:,i)),MED(:,i),pars.method) ],...
            [MED(:,i) MED(:,i)],medc,'LineWidth',2);
        catch
            medError = true;
        end
    end
 end
 %-------------------------------------------------------------------------
if pars.plotMean
    if ~medError
     L=legend([p(1) p(2)],'Mean','Median');
    else
        L=legend(p(1),'Mean');
    end
    set(L,'box','off','FontSize',14)
end

set(gca,'XTick',1:nData,'XTickLabels',xL);
 
 %-------------------------------------------------------------------------

if ~isempty(pars.medianStyle)
    meds = cellfun(@nanmedian,X);
    hold on;
    plot(1:length(X),meds,pars.medianStyle{:});
end


end %of function