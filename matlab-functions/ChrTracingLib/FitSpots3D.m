function [xyz_out,fitTable] = FitSpots3D(imIn,xyz_in,varargin)
% inputs
%    imIn - a 3D image
%    xyz_in - a list of 3D centroids, to pixel accuracy
% outputs
%    xyz - a list of 3D centroids to sub-pixel accuracy, Gaussian fit
%    fitTable - fit stats from FitPsf3D 

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'r_xy','integer',6};%  box radius in xy
defaults(end+1,:) = {'r_z','integer',4}; % box radius in z 
defaults(end+1,:) = {'showPlot','boolean',false};
defaults(end+1,:) = {'removeNonFits','boolean',true};

pars = ParseVariableArguments(varargin, defaults, mfilename);

[imBoxes,skipPts,xyz_in2] = BoxSpots3D(imIn,xyz_in,'parameters',pars);
nPts2 = length(imBoxes);
fitTable = cell(nPts2,1);
xyz_out =cell(nPts2,1); %  nan(nPts2,3);  % should parfor this
parfor b=1:nPts2
    dfit.x = [];
    if ~isempty(imBoxes{b})
        dfit = FitPsf3D(imBoxes{b},'seedPoint',[pars.r_xy+1,pars.r_xy+1,pars.r_z+1],'maxFitWidth',2*pars.r_xy+1); % 'troubleshoot',true,
        fitTable{b} = dfit;
    else
        x = []; y = []; z=[];
        dfit = table(x,y,z);
        xyz_out{b} = [nan,nan,nan];
    end
   if ~isempty(dfit.x)
       % xyz_out(b,:) = [xyz_in(b,1)-pars.r_xy-1+dfit.x, xyz_in(b,2)-pars.r_xy-1+dfit.y, xyz_in(b,3)-pars.r_z-1+dfit.z];
       xyz_out{b} = [xyz_in2(b,1)-pars.r_xy-1+dfit.x, xyz_in2(b,2)-pars.r_xy-1+dfit.y, xyz_in2(b,3)-pars.r_z-1+dfit.z];
   else
       xyz_out{b} = [nan,nan,nan];
   end
end
fitTable = cat(1,fitTable{:});
xyz_out = cat(1,xyz_out{:});
if pars.removeNonFits
    drop = isnan(xyz_out(:,1));
    xyz_out(drop,:) = [];
end


if pars.showPlot || nargout==0
    bMax = max(fitTable.a + fitTable.b);
    ProjectIm3D(imIn,'showPlots',true,'caxis',[0,bMax]);
    subplot(1,3,1); hold on; plot(xyz_in(:,1),xyz_in(:,2),'ro'); title('xy')
    subplot(1,3,2); hold on; plot(xyz_in(:,2),xyz_in(:,3),'ro'); title('yz')
    subplot(1,3,3); hold on; plot(xyz_in(:,1),xyz_in(:,3),'ro'); title('xz')
    subplot(1,3,1); hold on; plot(xyz_out(:,1),xyz_out(:,2),'r+'); title('xy')
    subplot(1,3,2); hold on; plot(xyz_out(:,2),xyz_out(:,3),'r+'); title('yz')
    subplot(1,3,3); hold on; plot(xyz_out(:,1),xyz_out(:,3),'r+'); title('xz')
end