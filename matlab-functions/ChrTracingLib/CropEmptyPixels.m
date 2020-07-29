function [im,boxCut] = CropEmptyPixels(im)

hasDataX = sum(im,1) > 0;
hasDataY = sum(im,2) > 1;
boxCut.xStart = find(hasDataX,1,'first');
boxCut.xStop = find(hasDataX,1,'last');
boxCut.yStart = find(hasDataY,1,'first');
boxCut.yStop = find(hasDataY,1,'last');
im = im(boxCut.yStart:boxCut.yStop,boxCut.xStart:boxCut.xStop,:);