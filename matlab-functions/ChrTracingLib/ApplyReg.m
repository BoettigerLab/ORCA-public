function imOut = ApplyReg(imIn,regData,varargin)
% imOut = ApplyReg(imIn,regData)
% 
%----------------
% Updates
%----------------
% 2018-11-24, updated to allow output from CorrAlignFast, which performs
% rotation before translation and uses theta for the angle instead of
% angle. We exploit this difference in the structure to preserve backwards
% compatibility.
%
% 2019-10-23  Updated to fix bug in rescaled images.
% CorrAlignFast performs scaling first, then rotation then
% translation. We need to do this in the same order. 
%
% 2021-08-17 Updated to fix bug resulting from 2-step save format for
% CorrAlignFast (2019-10-23), need to combine xshift2 and yshift2
% 
% see also:
%  CorrAlignRotateScale - compute a registration
%  CorrAlignFast - 2D registration performed at 2 scales for speed
%  RotateTranslatePoints - 2D rotation around a fixed point. optionally
%  rotate and translate an imagetiles cell array as well.
%  ImageTranslate - an accelerated version of Matlab's builtin imtranslate.
%                   (should speed test against Matlab 2020b). 

defaults = cell(0,3);
defaults(end+1,:) = {'rescaleShifts','float',1};
defaults(end+1,:) = {'invert','boolean',false};
defaults(end+1,:) = {'rotateMethod',{'bilinear','nearest'},'bilinear'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if istable(regData)
    regData = table2struct(regData);
end

% handle shift rescaling (mostly obsolete)
if isfield(regData,'xshift') && isfield(regData,'yshift') && pars.rescaleShifts ~= 1
    regData.xshift = regData.xshift*pars.rescaleShifts;
    regData.yshift = regData.yshift*pars.rescaleShifts;
    regData.rescale = 1-(1-regData.rescale)*pars.rescaleShifts;
end

% image rescaling is only 2D
if isfield(regData,'rescale')
    if regData.rescale ~= 1
        if size(imIn,3) > 1
            error('image rescaling only available for 2D images');
        end
    end
end

 % For speed in parsing output of CorrAlignFast
 %     notably, CorrAlignFast does rotation1, translation1, then rotation2
 %     translation2.  As long as rotation2 is zero, we can combine these
 %     translations for speed.  
    singleStep = true;
    if isfield(regData,'xshift2') && isfield(regData,'yshift2') 
        if regData.theta2 == 0 % if there's no rotation, we can combine these two translation steps
            singleStep = true;
            regData.xshift = regData.xshift + regData.xshift2;
            regData.yshift = regData.yshift + regData.yshift2;
        elseif regData.theta ~= 0
            singleStep = false; 
        end
    end
    
    
    imOut = imIn; 
    try
        if ~pars.invert
            % the original direction
            if isfield(regData,'angle') 
                % fxns that use "angle" instead of "theta" do translation
                % before rotation.  
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
                if regData.angle ~= 0 
                    imOut = imrotate(imOut, regData.angle,pars.rotateMethod,'crop');
                end
            elseif isfield(regData,'theta') % CorrAlignFast does (scale) then rot then trans
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                if regData.theta ~=0
                    imOut = imrotate(imOut, regData.theta,pars.rotateMethod,'crop');
                end
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
            else
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
            end
            % if 2 step with rotation, the second shift is peformed after 
            if ~singleStep
                twoStp.rescale = regData.rescale2;
                imOut = RescaleIm(imOut,twoStp,'parameters',pars);
                imOut = imrotate(imOut, regData.theta2,pars.rotateMethod,'crop');
                imOut = ImageTranslate(imOut,[regData.xshift2,regData.yshift2]);
            end
            
            
        else
            % the inverse direction
            if isfield(regData,'angle') % note, order also switches in inverse
                if regData.angle ~=0 
                    imOut = imrotate(imOut, -regData.angle,pars.rotateMethod,'crop');
                end
                imOut = ImageTranslate(imOut,-[regData.xshift,regData.yshift]);
                imOut = RescaleIm(imOut,regData,'parameters',pars);
            elseif isfield(regData,'theta') % note, order also switches in inverse
                imOut = ImageTranslate(imOut,-[regData.xshift,regData.yshift]);
                if regData.theta ~= 0
                    imOut = imrotate(imOut, -regData.theta,pars.rotateMethod,'crop');
                end
                imOut = RescaleIm(imOut,regData,'parameters',pars);
            else
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = ImageTranslate(imOut,-[regData.xshift,regData.yshift]);
            end
        end
        
        % apply imwarp if requested
        if isfield(regData,'tform')
            if ~isempty(regData.tform)
                imOut = imwarp(imOut,regData.tform.tobj,'OutputView',imref2d(size(imOut)));
            end
        end
        
    catch er
        warning(er.message);
        disp('place debug here');
        disp(er.getReport);
    end
end

function imOut = RescaleIm(imOut,regData,varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'invert','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);
% rescale matrix
    if isfield(regData,'rescale')
        if regData.rescale ~= 1 
            sc =regData.rescale;
            if pars.invert
                sc=1/sc;
            end
            S = [sc   0   0;
                 0   sc  0;
                 0   0   1];
            imOut = imwarp(imOut,affine2d(S) ,'OutputView',imref2d(size(imOut))  );
        end
    end
end