function daxNames = SaveAlignedDax(fidMapData,regData,varargin)
%% Apply registration and save new images

defaults = cell(0,3);
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'tag','string',''};
defaults(end+1,:) = {'fov','integer',0};
defaults(end+1,:) = {'dataType','string','big endian'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

numHybes = size(fidMapData,1);
daxNames = cell(numHybes,1);

for h=1:numHybes
    xshift = regData(h).xshift;
    yshift = regData(h).yshift;
    angle = regData(h).angle;
    tform = regData(h).tform;
    % rotation must be applied to whole image, not to part. 
    temp = ReadFromMemMap(fidMapData{h}); %  figure(10); clf; imagesc(max(temp,[],3));
    temp = imrotate(ImageTranslate(temp,[xshift,yshift]),angle,'bilinear','crop');  % imtranslate
    if ~isempty(tform)
        temp = imwarp(temp,tform.tobj,'OutputView',imref2d(size(temp)));
    end
    name = ['fov',num2str(pars.fov,'%03d'),'_h',num2str(h,'%03d'),'_',pars.tag];
    WriteDax(temp,'daxName',name,'folder',pars.saveFolder,'dataType',pars.dataType);
    daxNames{h} = [name,'.dax'];
end   