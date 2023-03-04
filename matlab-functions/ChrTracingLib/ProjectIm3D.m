function [xy,xz,yz] = ProjectIm3D(im3D,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'caxis','array',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

xy = max(im3D,[],3);
yz = max( permute(im3D,[3,1,2]), [],3);
xz = max( permute(im3D,[3,2,1]), [],3);

if pars.showPlots || nargout==0
   subplot(1,3,1); imagesc(xy); colorbar; if ~isempty(pars.caxis); caxis(pars.caxis); end;
   subplot(1,3,2); imagesc(yz); colorbar;  if ~isempty(pars.caxis); caxis(pars.caxis); end;
   subplot(1,3,3); imagesc(xz); colorbar;  if ~isempty(pars.caxis); caxis(pars.caxis); end;
end