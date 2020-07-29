function cellData = PlotDaxImages(imageData,varargin)

%%

% Get image dimensions
daxfile = regexprep(imageData(1).name,['_',imageData(1).binType,'\.bin'],'\.dax');
fileFolder = [fileparts(imageData(1).filePath),filesep];
infFile = ReadInfoFile([fileFolder,daxfile]); 
H = infFile.frame_dimensions(1); 
W = infFile.frame_dimensions(2); 

% Get number of hybes
numHybes = length([imageData.hybNum]);


% -------------------------------------------------------------------------
% Plotting transformed data
% -------------------------------------------------------------------------
alignedDaxID = find(strcmp('alignedDaxImages', parameters.reportsToGenerate(:,1)));
if ~isempty(alignedDaxID)
    alignedDaxHandle = figure('Name',['alignedDaxImages_cell', num2str(imageData(1).cellNum)], ...
        'visible', parameters.reportsToGenerate{alignedDaxID, 2});
end

tiledDaxID = find(strcmp('tiledDaxImages', parameters.reportsToGenerate(:,1)));
if ~isempty(tiledDaxID)
    tiledDaxHandle = figure('Name',['tiledDaxImages_cell', num2str(imageData(1).cellNum)], ...
        'visible', parameters.reportsToGenerate{tiledDaxID, 2});
    set(gcf,'Units','Inches','Position',[1,1,2.36*numHybes/4,9.5],'PaperUnits', 'Inches', 'PaperSize', [2.36*numHybes/4, 9.5]);
end


if ~isempty(tiledDaxID) || ~isempty(alignedDaxID)
    rawIm  = zeros(H,W,numHybes);
    alignedIm = zeros(H,W,imageData.cellNum,'uint16');
    for h=1:numHybes
        daxfile = regexprep(imageData(h).name,['_',imageData(h).binType,'\.bin'],'\.dax');
        dax = max(ReadDax([fileFolder,daxfile],'endFrame',1),[],3); % max project the first frame
        tformInv = fliptform(imageData(h).tform); 
        alignedDax = imtransform(dax,tformInv,...
                        'XYScale',1,'XData',[1 W],'YData',[1 H]);
        rawIm(:,:,h) = dax;
        alignedIm(:,:,h) = uint16(alignedDax);
        if ~isempty(tiledDaxID)
            figure(tiledDaxHandle); subplot(4,numHybes/4,h); 
            imagesc(uint16(alignedDax)); colormap gray;  % caxis([0,2^13.5]);
        end
    end
end

% save alignedDax image if requested
if ~isempty(alignedDaxID)
    figure(alignedDaxHandle);
    Ncolor(alignedDax,jet(numHybes));
    FluorImage('box',false);
    SaveFigure(alignedDaxHandle,'overwrite',parameters.overwrite,'formats','.png');
    close(alignedDaxHandle); 
end

% save tileDax image if requested
if ~isempty(tiledDaxID)   
    figure(tiledDaxHandle)
    FluorImage('box',false);
    SaveFigure(tiledDaxHandle,'overwrite',parameters.overwrite,'formats',parameters.figFormats);
    close(tiledDaxHandle); 
end


% export data
cellData.alignedDax = alignedIm; 
cellData.imageSize = [H,W];
cellData.numHybes = numHybes;

