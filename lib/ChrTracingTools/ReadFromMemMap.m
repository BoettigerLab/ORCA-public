function movie = ReadFromMemMap(memMapData,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'startFrame','integer',1}; % 0 = first frame
defaults(end+1,:) = {'endFrame','integer',[]}; % 0 = last frames
defaults(end+1,:) = {'selectFrames','integer',[]}; % 0 = all frames
defaults(end+1,:) = {'roi','integer',[]}; % 0 = full field
defaults(end+1,:) = {'mapName','string',''};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% parse 
framesInDax = memMapData.framesInDax;
frameDim = memMapData.frameDim;
if isempty(pars.selectFrames)
    selectFrames = memMapData.selectFrames;
else
    selectFrames = pars.selectFrames;
end
framesToLoad = length(selectFrames);

if isempty(pars.roi)
   xi = 1; xe = frameDim(1);
   yi = 1; ye = frameDim(2);
else
    xi = pars.roi(1);
    xe = pars.roi(2);
    yi = pars.roi(3);
    ye = pars.roi(4);
end

[ri,ci,zi] = meshgrid(uint32(xi:xe),uint32(yi:ye),selectFrames);
inds = sub2indFast([frameDim(2),frameDim(1),framesInDax],ri(:),ci(:),zi(:));
movie = memMapData.memMap.Data(inds);
if strcmp(memMapData.binaryFormat,'b')
    movie = swapbytes(movie);
end
xs = xe-xi+uint32(1);
ys = ye-yi+uint32(1);
movie = reshape(movie,[ys,xs,framesToLoad]); % figure(7); clf; imagesc(max(movie,[],3)); colormap gray;

