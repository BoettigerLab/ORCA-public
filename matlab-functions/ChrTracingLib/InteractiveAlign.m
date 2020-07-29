function [alignValues,comboImage,valid,pars] = InteractiveAlign(refImage,addImage,varargin)
% Merge "addImage" into existing "refImage" using automated (CorrAlignFast)
%   and interactive correction.
% 
% 
%%
global MVR

% ------------ Parse default Parameters
defaults = cell(0,3);
% InteractiveAlign pars
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'startPosition','array',[0,0]};
defaults(end+1,:) = {'combine',{'overlay','merge','off'},'overlay'};
defaults(end+1,:) = {'outputTform','boolean',true}; % make default or not  
defaults(end+1,:) = {'refImageCnst','boolean',true}; % accelerate by passing an already contrasted ref Image;  
defaults(end+1,:) = {'imName','string','1'}; 
% MosaicViewerRender pars
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur'},'edgeBlur'};
defaults(end+1,:) = {'downsample','integer',1};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'interactive','boolean',true};
defaults(end+1,:) = {'minOverlap','fraction',.05};
defaults(end+1,:) = {'corrFig','integer',9}; % figure to use for correlation plot 
defaults(end+1,:) = {'interactFig','integer',30}; % figure to use for interactive plot 
defaults(end+1,:) = {'xyShifts','freeType',[]}; %
defaults(end+1,:) = {'padMosaic','nonnegative',3}; % in multiples of tile size
% parameters for prepping CorrAlign
defaults(end+1,:) = {'corrAlign','boolean',true};
defaults(end+1,:) = {'padAlign','integer',150}; % pixels to add on either side of reference image
defaults(end+1,:) = {'corrAlignHigh','boolean',.9999};
defaults(end+1,:) = {'corrAlignLow','boolean',.8};
defaults(end+1,:) = {'sortData','boolean',true};
% parameters forwarded to CorrAlignFast
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % 0.9:0.01:1.10
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'saveCorr', 'boolean', false};
defaults(end+1,:) = {'maxSize', 'positive', 300};
defaults(end+1,:) = {'maxShift', 'nonnegative', inf}; % 0= auto to 1/5th image height 
defaults(end+1,:) = {'minGrad', 'float', -inf}; % 50
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'savePath', 'string', 'G:\Alistair\CorrAlign\'};
defaults(end+1,:) = {'showExtraPlot','boolean',false};
defaults(end+1,:) = {'fineUpsample', 'positive', 1}; 
defaults(end+1,:) = {'fineAngles','float',0}; % -10:1:10
defaults(end+1,:) = {'fineScales','float',1}; % 0.9:0.01:1.10
defaults(end+1,:) = {'fineCenter','array',[0,0]}; % 0 0 = use brightest
defaults(end+1,:) = {'fineBox', 'freeType', []}; 
defaults(end+1,:) = {'minFineImprovement', 'float', .1};
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);

if pars.interactive || ~strcmp(pars.combine,'off') 
    figH = figure(pars.interactFig);   pause(.01);
    MVR.figH = figH;
end
if pars.interactive
    set(figH,'KeyPressFcn',{@FigKeyPress});
    MVR.im3 = zeros(10); 
    ResetShifts();
end


% initialize a empty alignValues
alignValues.xshift = 0; 
alignValues.yshift = 0;
alignValues.theta = 0;
alignValues.rescale = 1;

[h_i,w_i] = size(addImage);
if pars.maxShift == 0 
   pars.maxShift = h_i/5;
end
if pars.maxShift < 1
    pars.maxShift = h_i*pars.maxShift;
end

uls = pars.startPosition;
a = pars.padAlign;
m=1; % 1 image at a time, fixed for now

%--- Stp 1: crop a region from the reference image. 
% protect against edge out-of-bounds errors
boxCoords(m,:) = round([uls(m,1),uls(m,1)+w_i-1,uls(m,2),uls(m,2)+h_i-1]);
[h_m,w_m] = size(refImage);
% don't even try to place tiles that are outside of the box
noOverlap = false;
if boxCoords(m,3) < -h_i+2*a ...
        || boxCoords(m,4) > h_m+h_i-2*a ...
        || boxCoords(m,1) < -w_i+2*a ...
        || boxCoords(m,2) > w_m+w_i-2*a
    noOverlap = true;
end
if noOverlap
    blnkImage = zeros(size(addImage),class(addImage));
    comboImage = cat(3,blnkImage,addImage);
    valid = 4; % no use coming back to this one;
    return
end
    
y1 = max(boxCoords(m,3)-a,1);
x1 = max(boxCoords(m,1)-a,1);
y2 = min(boxCoords(m,4)+a,h_m);
x2 = min(boxCoords(m,2)+a,w_m);
currImage = refImage(y1:y2,x1:x2);
xOff = 0;
yOff = 0;
% protect against image size difference errors induced by edge out-of-bounds 
if x1==1
    xOff = -(boxCoords(m,1)-a)+1;
    currImage = padarray(currImage,[0,xOff],nan,'pre');
end
if x2==w_m
    currImage = padarray(currImage,[0,(boxCoords(m,2)+a)-w_m],nan,'post');
end
if y1==1
    yOff = -(boxCoords(m,3)-a)+1;
   currImage = padarray(currImage,[yOff,0],nan,'pre'); 
end
if y2==h_m
    currImage = padarray(currImage,[(boxCoords(m,4)+a)-h_m,0],nan,'post');
end

% adjust contrast
fractionHasData = nansum(currImage(:)>0)./length(currImage(:));
low = (1-fractionHasData)+pars.corrAlignLow*fractionHasData; 
high = 1-(1-pars.corrAlignHigh)*fractionHasData;
MVR.contrasts.refMax = high;
MVR.contrasts.refMin = low;
if pars.refImageCnst
    currImageCnst = IncreaseContrast(currImage,'low',low,'high',high);
else
    currImageCnst = currImage;
end
addImageCnst = IncreaseContrast(addImage,'low',pars.corrAlignLow,'high',pars.corrAlignHigh);
addImageCnst = padarray(addImageCnst,[a,a],nan);
MVR.contrasts.datMax = pars.corrAlignHigh;
MVR.contrasts.datMin = pars.corrAlignLow;

% figure(3); clf; 
% subplot(2,2,1); imagesc(refImage);
% subplot(2,2,2); imagesc(addImage);
% subplot(2,2,3); imagesc(currImageCnst);
% subplot(2,2,4); imagesc(addImageCnst);
% figure(4); clf; Ncolor(cat(3,currImageCnst,addImageCnst));





valid = 0; 
% 0 - image not aligned
% 1 - image aligned
% 2 - skip image, to align later
% 3 - flag image (for whatever)
% 4 - no Ref image data - use a different alignment (from file)

% if reference image is blank, leave the add image where it was
if sum(currImageCnst(:)>0) < 10 % or almost blank
    pars.corrAlign = false;
    pars.interactive=false;
    valid = 4;
end



if pars.corrAlign
   % use corrAlign to place tile
   if fractionHasData > pars.minOverlap %#ok<*BDSCI>            
        if pars.showplot 
            figOut = figure(pars.corrFig); clf; 
        end
        alignValues = CorrAlignFast(currImageCnst,addImageCnst,...
            'parameters',pars); % 
        if pars.showplot
            subplot(2,2,1); title(pars.imName);
            SaveFigure(figOut,'name',['corrFig_',pars.imName],...
                'formats',{'png'},'overwrite',true,'saveData',pars.saveCorr);
        end
   end
end

if pars.interactive && fractionHasData > pars.minOverlap 
    figure(100); clf; 
    Ncolor(cat(3,currImageCnst,addImageCnst));    
    title(['Original image ', pars.imName]);  
    addImageCnst_align = ScaleRotateShift(addImageCnst,alignValues); % 
    im3 = cat(3,currImageCnst,addImageCnst_align);
    figure(MVR.figH); clf; Ncolor(im3); 
     title(['Aligned image ', pars.imName,'. adjust manually if necessary ']);  
    if pars.verbose
        disp('auto computed values');
        disp(alignValues);  % added
        % disp(MVR.shifts);
    end
    MVR.im3 = im3;
    textDir = {'Use keyboard to adjust alignment:';
        '"a"/"d" for left/right, "w"/"s" for up/down';
        '"j"/"l" for left/right fast,"i/k" for up down fast';
        '"e"/"r" to rotate cw/ccw, "o"/"p" to rotate fast';
        '"y"/"u" to expand/contract';
        '"v"/"b" increase/decrease red contrast';
        '"n"/"m" increase/decrease cyan contrast';
        'To manually select control points instead, press "c"';
        'To cancel, press "x", to skip for now, press "z"';
        'Click "Okay" or close this box to save and continue.'};
    okayPressed = msgbox(textDir);
    set(okayPressed,'KeyPressFcn',{@FigKeyPress});
    set(MVR.figH,'KeyPressFcn',{@FigKeyPress});          
    waitfor(okayPressed);
    if ~MVR.shifts.skip && ~MVR.shifts.cancel && ~MVR.shifts.useControlPoints
        disp(MVR.shifts);
        alignValues.xshift = alignValues.xshift + MVR.shifts.leftRightOffset;
        alignValues.yshift = alignValues.yshift + MVR.shifts.upDownOffset;
        alignValues.theta = alignValues.theta + MVR.shifts.angle;
        alignValues.rescale = alignValues.rescale*MVR.shifts.rescale;
        valid = 1;
    elseif MVR.shifts.useControlPoints
        [selectedMovingPoints,selectedFixedPoints] = cpselect(addImageCnst,currImageCnst,'Wait',true);
        tformCP = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'affine');
        warp = WarpProps(tformCP,'imsize',size(addImageCnst)); % not entirely sure this is accurate. 
        alignValues.xshift = alignValues.xshift + warp.xshift;
        alignValues.yshift = alignValues.yshift + warp.yshift;
        alignValues.theta = alignValues.theta + warp.rotationAngle;
        alignValues.rescale = alignValues.rescale*warp.scaleImage;
        valid =1;
        if pars.verbose
            [cx,cy]=transformPointsInverse(tformCP,selectedMovingPoints(:,1),selectedMovingPoints(:,2));
            errWarp = mean(sqrt(sum(([cx,cy] - selectedFixedPoints).^2,2)));
            disp(['mean warp error =',num2str(errWarp),' pixels']);
        end
    elseif MVR.shifts.skip
        valid = 2; % delay
    elseif MVR.shifts.cancel
        valid = 0;
    end
    if MVR.shifts.flag && ~MVR.shifts.skip
        valid = 3; % keep but mark
    end
    ResetShifts();
else
    valid = 4; % skip this completely if insufficient overlap for alginment
    alignValues.xshift = 0; % and keep alignValues at 0
    alignValues.yshift = 0;
    alignValues.theta = 0;
    alignValues.rescale = 1;
end

if ~strcmp(pars.combine,'off')
    try
    % just for viewing. The main goal is to compute the shifts. 
    comboImage =currImageCnst; % refImage;
%     addImagePad = padarray(addImage,[a,a],nan);
%     alignValues
    addImage2 = ScaleRotateShift(addImageCnst,alignValues); % switched to ApplyReg from ScaleRotateShift 10/23/19 
    if strcmp(pars.combine,'overlay')
        comboImage = cat(3,currImageCnst,addImage2);
    elseif strcmp(pars.combine,'merge')
        comboImage= CombineImages(comboImage,addImage2,'method',pars.method);
    end
    figure(MVR.figH); clf; Ncolor(IncreaseContrast(comboImage,'high',.9995));
    if pars.verbose
        disp('saved values'); disp(alignValues);
    end
    % figure(31); clf; Ncolor(IncreaseContrast(MVR.im3,'high',.9995)); title('aligned');
    % figure(32); clf; Ncolor(IncreaseContrast(comboImage,'high',.9995)); title('saved');
    % figure(33); clf;  Ncolor(IncreaseContrast(cat(3,currImageCnst,addImage2),'high',.9995)); title('orig');
    catch er
        disp(er.getReport);
        error('something wrong in InteractiveAlign Combine images');
    end
else
    comboImage = [];
end


%%
% --- called when figureH is in view and a key is pressed
function FigKeyPress(hObject,eventdata,handles)  %#ok<INUSD,INUSL>
    global MVR
    key = eventdata.Key;
    shifts = MVR.shifts; % could also have stored this in handles
    contrasts = MVR.contrasts;
    % can we be less lazy and make this a pop-up GUI to compliment the long
    % list of keyboard commands?
    switch(key)
        case('w') %want to move top image up
            shifts.upDownOffset = shifts.upDownOffset - 1;
        case('s') %want to move top image down
            shifts.upDownOffset = shifts.upDownOffset + 1;
        case('a') %want to move top image left
            shifts.leftRightOffset = shifts.leftRightOffset - 1;
        case('d') %want to move top image right
            shifts.leftRightOffset = shifts.leftRightOffset + 1;
        case('i') %want to move top image up fast
            shifts.upDownOffset = shifts.upDownOffset - 10;
        case('k') %want to move top image down fast
            shifts.upDownOffset = shifts.upDownOffset + 10;
        case('j') %want to move top image left fast
            shifts.leftRightOffset = shifts.leftRightOffset - 10;
        case('l') %want to move top image right fast
            shifts.leftRightOffset = shifts.leftRightOffset + 10;
        case('e') % rotate small step CW
            shifts.angle = shifts.angle + .25;
        case('r') % rotate small step CCW
            shifts.angle = shifts.angle - .25;
        case('o') % rotate large step CW
            shifts.angle = shifts.angle + 2.5;
        case('p')% rotate large step CCW
            shifts.angle = shifts.angle - 2.5;
        case('y') % expand
            shifts.rescale = shifts.rescale + .001;
        case('u') % contract
            shifts.rescale = shifts.rescale - .001;
        case('v') % max contrast ref down
            a = contrasts.refMax;
            contrasts.refMax = a-(1-a)*.1;
        case('b') % max contrast ref up
            a = contrasts.refMax;
            contrasts.refMax = a+(1-a)*.1;
        case('n') % max contrast data down
            a = contrasts.datMax;
            contrasts.datMax = a-(1-a)*.1;
        case('m') % max contrast data up
            a = contrasts.datMax;
            contrasts.datMax = a+(1-a)*.1;
        case('c')
            shifts.useControlPoints = true;
        case('x')
            shifts.cancel = true;
        case('z')
            shifts.skip = true;
        case('f')
            shifts.flag = true;
        case(28) %For right and left arrows,
        case(29) %will add this soon
            % should add an image 
    otherwise
    end
    MVR.shifts = shifts; % could also stick into handles;
    MVR.lastKey = key;
    MVR.contrasts = contrasts;
    UpdateOverlay();

function UpdateOverlay
    global MVR
    % can we get figure zoom and set figure zoom (not just a plot)
    % can we add the im1 im2 high/low to the keypad?
    im1 = MVR.im3(:,:,1);
    im2 = MVR.im3(:,:,2);
    aVs.xshift = MVR.shifts.leftRightOffset;
    aVs.yshift = MVR.shifts.upDownOffset;
    aVs.theta = MVR.shifts.angle;
    aVs.rescale = MVR.shifts.rescale;
    im2 = ScaleRotateShift(im2,aVs); % swapped
    im1 = IncreaseContrast(im1,'low',MVR.contrasts.refMin,'high',MVR.contrasts.refMax);
    im2 = IncreaseContrast(im2,'low',MVR.contrasts.datMin,'high',MVR.contrasts.datMax);
    im3 = cat(3,im1,im2);
    figure(MVR.figH); cla;
    Ncolor(im3); 

function ResetShifts
    global MVR
    MVR.shifts.leftRightOffset = 0;
    MVR.shifts.upDownOffset = 0;
    MVR.shifts.angle = 0;
    MVR.shifts.rescale = 1;
    MVR.shifts.skip = false;
    MVR.shifts.cancel = false;
    MVR.shifts.flag = false;
    MVR.shifts.useControlPoints = false;