%% a good object to fit Z on?
function currStack = FitZ(currStack,varargin)
%  examines the brightness in each step of the z-stack and fits these data
%  to a precomputed gaussian PSF to determine the center. 
%  If there is not enough data to compute or if the fit quality is bad, the
%  function will return NaN for the position. 
% 
% 
%% Inputs  
% currStack - a 3D matrix, t_Obs x d_Dims x z_Depth.   
%             d-Dims are x,y,z,t,h,b: x-position, y-position, z-position,
%             t frame number, gaussian-fit height and gaussian fit
%             background.  
%             t-Obs is number of total z-stacks (not total frames)
%             
%% Outputs
% currStack - updated version of the input matrix
%           - the z-position data (typically nan) will be filled in with nm
%             computed z-position based on the brightness of the spot in
%             the recorded positions in the z-stack. 

%% parse optional parameters
defaults = cell(0,3);

defaults(end+1,:) = {'showFig','integer',10};
defaults(end+1,:) = {'ciMax','fraction',0.75}; % if this confidence interval range is larger than z-range, the fit will be rejected.  0.95 is stringent, 0.75 less stringent, 0.5 relaxed
defaults(end+1,:) = {'spotNum','integer',nan};
defaults(end+1,:) = {'displayInterval','boolean',500};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);



% more parameters 
[tObs,dDims,zDepth] = size(currStack); 
% dDims:  x,y,z,t,h,b

z = (1000:1000:(1000*(zDepth)))'; % hard-coded 1000 nm steps
%  This goes hand-in-hand with hard-coded z-PSF height

ftype = fittype('exp(-((x-mu_x)/(2*500)).^2)*a',...
                 'coeff', {'a','mu_x'},'ind',{'x'});
% simple Gaussian with 0 background (since spot heights are already
% background subtracted).  The width of the PSF in Z is fixed (500 nm)



%% start fitting
k=0;
currStack(:,3,:) = nan; 

for t=1:tObs
    b = squeeze(currStack(t,5,:)); % column 5 is "brightness"
    if  any(diff(b)) % at least 2 consecutive data points needed for a fit
        b(isnan(b)) = 0;
        % disp('fitting')
         try
           %    Gaussian fit
           [~,i] = max(b);
           sp = [max(b),    z(i)];
           lb = [   0,         0]; 
           ub = [4*max(b),  max(z)];
    
            fitZ = fit(z,b,ftype,...
                 'StartPoint',sp,...
                 'Lower',lb,...
                 'Upper',ub,...
                 'TolFun',1e-6,'TolX',1e-6,'MaxIter',4000,'MaxFunEvals',600);
           %  figure(10); clf; plot(fitZ,z,b); pause(.1); xlabel('z'); ylabel('brightness'); 
            currStack(t,3,:) = repmat(fitZ.mu_x,zDepth,1);
            k=k+1;

            ci = confint(fitZ,.75);
            if (ci(2,2)-ci(1,2))/(max(z)-min(z)) > 1
                poorSpot = ['detected poor fit on spot = ',num2str(pars.spotNum),' t=',num2str(t)];
                if pars.veryverbose
                    disp(poorSpot);
                end
                
                currStack(t,3,:) = nan;
                % figure(pars.showFig); clf; 
                % plot(fitZ,z,b); pause(.1); 
                % xlabel('z'); ylabel('brightness'); title(poorSpot);
                % disp('place debug here')
            end
            if  rem(k,pars.displayInterval) == 0 && pars.showFig
              figure(pars.showFig); clf; plot(fitZ,z,b); 
              pause(.1); xlabel('z (nm)'); ylabel('brightness'); 
              title(['example fit on spot = ',num2str(pars.spotNum),' t=',num2str(t)]);
              pause(.01); 
            end

        catch er
            disp(er.getReport);
            disp('place debug here')
         end  
    end
end

if pars.showFig
    xs = 1:zDepth:(zDepth*tObs);
    ys = squeeze(currStack(:,3,1));
    figure(pars.showFig); clf; 
    plot(xs(~isnan(ys)),ys(~isnan(ys)),'.-');
     % 
     % if length(ys(~isnan(ys))) > 300
     %     disp('test')
     % end

    ylim([.8*min(z),1.1*max(z)])
    xlim([1,(zDepth*tObs)]); 
    xlabel('time (frame)');
    ylabel('spot z (nm)');
    aveZstep =  median(abs(diff(ys(~isnan(ys))))) ;
    title(['median z-step per observation:'  ,num2str(aveZstep,3), ' nm ' ]);
    aveZstep = nanmedian(abs(diff(ys))) ;
    title(['Spot ',num2str(pars.spotNum), ' median z-step per ', num2str(zDepth), ' frames: '  ,num2str(aveZstep,3), ' nm ' ]);   
end
