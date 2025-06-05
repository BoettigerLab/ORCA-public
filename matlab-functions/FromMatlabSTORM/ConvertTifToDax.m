function ConvertTifToDax(tifFiles)
% Convert Hal large Tif to Dax. 
% 
% NAS02_Vol4 = 'L:\';
% dataFolder = [NAS02_Vol4,'Jude\2023-11-30_50kb\'];
% tifFiles = FindFiles([dataFolder,'beads*.tif']); blockSize = 99; nFrames = 99;
% tifFiles = FindFiles([dataFolder,'sample*.tif']); blockSize = 100; nFrames = 18000;
% 
% Alistair Boettiger, April 2024
defaultBlockSize = 100;
nW = length(tifFiles);
for w=1:nW 
    tic
    disp(['coverting ',num2str(w),' of ',num2str(nW),' files']);
    xmlName =  regexprep(tifFiles{w},'.tif','.xml');
    xmlName = regexprep(xmlName,'_C1',''); % xml files don't contain camera labels, they recording setting for both cameras.
    xmlFile = ReadXML(xmlName);
    nFrames = xmlFile.settings.film.frames;
    blockSize = min(nFrames,defaultBlockSize);

    daxName = regexprep(tifFiles{w},'.tif','.dax');
    if exist(daxName,'file')
        disp('existing dax file found, skipping');
    else
        % read in the first frame block and start writing the dax
        im = tiffreadVolume(tifFiles{w},'PixelRegion',{[1 inf],[1 inf],[1 blockSize]});
        WriteDax(im,'saveFullName',daxName);
        infoFi = ReadInfoFile(daxName);
        infoFi.number_of_frames = nFrames;
        WriteInfoFiles(infoFi);
        pause(.1);
        % read in the subsequent frame blocks
        for k = 2:(nFrames/blockSize)
            firstFrame = (k-1)*blockSize+1;
            lastFrame = min(nFrames,k*blockSize);
            im = tiffreadVolume(tifFiles{w},'PixelRegion',...
                {[1 inf],[1 inf],[firstFrame, lastFrame]});
            AppendToDax(daxName,im);
            pause(.1);
        end
        tfin = toc;
        disp(['conversion complete in ',num2str(tfin/(60*60),3),' hours']);
    end
end

%% test view of new dax files
% daxS = ReadDax(daxName,'startFrame',1,'endFrame',2);
% daxE = ReadDax(daxName,'startFrame',nFrames-1,'endFrame',nFrames);
% figure(1); clf; 
% subplot(1,2,1); imagesc(IncreaseContrast( max(daxS,[],3),'high',.9999,'low',.3));
% subplot(1,2,2); imagesc(IncreaseContrast( max(daxE,[],3),'high',.9999,'low',.3));

