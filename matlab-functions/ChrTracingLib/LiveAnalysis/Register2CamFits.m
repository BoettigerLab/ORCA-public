function fitsOut = Register2CamFits(fitsIn,varargin)
% Apply camera registration to the fits structure (from LoadHD5Fits) 
% Apply chormatic correction to the registered image if required. 


defaults = cell(0,3);
% data loading
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'alignment_file','freeType',''}; % f
defaults(end+1,:) = {'chrom_correct_file','string',''}; % 
defaults(end+1,:) = {'nppXY','positive',108}; % 
pars = ParseVariableArguments(varargin,defaults,mfilename);
% matching

if ~isempty(pars.alignment_file) && (ischar(pars.alignment_file) || isstring(pars.alignment_file))
    alignT  = readtable(pars.alignment_file);
    alignS = table2struct(alignT);
    applyCameraAlign = true; % better readability
elseif ~isempty(pars.alignment_file) && isstruct(pars.alignment_file)
   alignS = pars.alignment_file;
    applyCameraAlign = true; % better readability
else
    alignS = [];
    applyCameraAlign = false;
end

if ~isempty(pars.chrom_correct_file)
    load(pars.chrom_correct_file)
    applyChromaticCorrect = true;
    if exist('tformPars','var')
        nppXY = tformPars.npp;
    else
        nppXY = pars.nppXY;
    end
else
    applyChromaticCorrect = false;
end

% === Apply camera alignment and chromatic correction
fitsOut = fitsIn;  % structures are a pretty good shape for this data
nFrames = length(fitsIn);
for f=1:nFrames
    if applyCameraAlign && strcmp(alignS.flip,'flipLR_none_none')
        pts2 = [alignS.w-fitsIn(f).x, fitsIn(f).y];
    else
        pts2 = [fitsIn(f).x, fitsIn(f).y];
    end
    if ~isempty(pts2)
        if applyCameraAlign
            pts2a = RotateTranslatePoints(pts2,alignS,'center',[alignS.h/2,alignS.w/2],'invert',false); 
        else
            pts2a = pts2;
            disp('skipping camera alignment')
        end
        if applyChromaticCorrect
            pts2b = pts2a*nppXY;
            pts2b = tforminv(tform3D,pts2b(:,1),pts2b(:,2),zeros(size(pts2b,1),1))/nppXY;  % not currently warping Z anyway
        else
            pts2b = pts2a;
        end
        fitsOut(f).x = pts2b(:,1);
        fitsOut(f).y = pts2b(:,2); 
    end
end
