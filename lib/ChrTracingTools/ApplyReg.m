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
    % 
    % see also:
    %  CorrAlignRotateScale - compute a registration
    %  CorrAlignFast
    %  RotateTranslatePoints - 2D rotation around a fixed point. optionally
    %  rotate and translate an imagetiles cell array as well.
    % 

    defaults = cell(0,3);
    defaults(end+1,:) = {'rescaleShifts','float',1};
    defaults(end+1,:) = {'invert','boolean',false};
    pars = ParseVariableArguments(varargin,defaults,mfilename);

    % 
    if isfield(regData,'xshift') && isfield(regData,'yshift') && pars.rescaleShifts ~= 1
        regData.xshift = regData.xshift*pars.rescaleShifts;
        regData.yshift = regData.yshift*pars.rescaleShifts;
        regData.rescale = 1-(1-regData.rescale)*pars.rescaleShifts;
    end


    imOut = imIn; 
    try
        if ~pars.invert
            % the original direction
            if isfield(regData,'angle') % fxns that save angle do trans then rot
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
                imOut = imrotate(imOut, regData.angle,'bilinear','crop');
            elseif isfield(regData,'theta') % CorrAlignFast does rot then trans
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = imrotate(imOut, regData.theta,'bilinear','crop');
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
            else
                imOut = RescaleIm(imOut,regData,'parameters',pars);
                imOut = ImageTranslate(imOut,[regData.xshift,regData.yshift]);
            end
        else
            % the inverse direction
            if isfield(regData,'angle') % note, order also switches in inverse
                imOut = imrotate(imOut, -regData.angle,'bilinear','crop');
                imOut = ImageTranslate(imOut,-[regData.xshift,regData.yshift]);
                imOut = RescaleIm(imOut,regData,'parameters',pars);
            elseif isfield(regData,'theta') % note, order also switches in inverse
                imOut = ImageTranslate(imOut,-[regData.xshift,regData.yshift]);
                imOut = imrotate(imOut, -regData.theta,'bilinear','crop');
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