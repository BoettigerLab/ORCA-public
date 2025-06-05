function im = ScaleRotateShift(im,alignValues,varargin)
% im = ScaleRotateShift(im,alignValues)
% 
%  OBSOLETE -- see ApplyReg
%% Inputs
% image im
% alignValues.rescale
% alignValues.theta
% alignValues.xshift
% alignValues.yshift
% 
%% Optional inputs
% 'opOrder', {'Scale','Rotate','Translate'} 
%    - change the order of operations using a cell of strings.  
% 'invert',true   
%    - invert the operation orginally specified by alignValues
%      This changes the direction of translation, inverts the rescaling,
%      and changes the direction of rotation. It does not invert the
%      'rescale shifts' value if passed.
% 'rescale shifts'
%    - rescale the scaling and translation by this value. This is helpful
%    if align values was computed on a downsampled image but you want to
%    apply the transformation to the full original image. 
% 
%% output
% returns the image after rescaling, rotating by theta, and shifting by
% xshift and yshift in that order. 
%
% Updates
% rewritten from scratch 2019-10-23, preserving backwards compatability
% new version allows order of operations to be changed more elegantly. 
% 
% see also ApplyReg

    defaults = cell(0,3);
    defaults(end+1,:) = {'rescaleShifts','float',1};
    defaults(end+1,:) = {'opOrder','cell',{'Scale','Rotate','Translate'}};
    defaults(end+1,:) = {'invert','boolean',false};
    pars = ParseVariableArguments(varargin,defaults,mfilename);

    % determine the order of operations requested
    opOrder = StringFind(pars.opOrder,{'Scale','Rotate','Translate'});
    if pars.invert % invert if requested
        opOrder = fliplr(opOrder); 
    end
    moveFxns = {@ImScale,@ImRot,@ImTrans};
    moveFxns = moveFxns(opOrder); 
    
    % apply the transforms (specified below) 
    if ~pars.invert
        im = moveFxns{3}(moveFxns{2}(moveFxns{1}(im,alignValues,'parameters',pars),alignValues,'parameters',pars),alignValues,'parameters',pars);
    end
    
    % ---sequential (non-commutive) transform
    % first we repopulate the align values
    if isfield(alignValues,'xshift2')
        av.xshift = alignValues.xshift2;
        av.yshift = alignValues.yshift2;
    end
    if isfield(alignValues,'theta2')
        av.theta = alignValues.theta2;
    end
    if isfield(alignValues,'rescale2')
        av.rescale = alignValues.rescale2;
    end
    
    if isfield(alignValues,'xshift2')
        im = moveFxns{3}(moveFxns{2}(moveFxns{1}(im,av,'parameters',pars),av,'parameters',pars),av,'parameters',pars);
    end
    
    % if inverting, we need to do the second set of shifts first (as well
    % as reverse the order of scale,rotate,shift)
    if pars.invert
        im = moveFxns{3}(moveFxns{2}(moveFxns{1}(im,alignValues,'parameters',pars),alignValues,'parameters',pars),alignValues,'parameters',pars);
    end
end

function im = ImScale(im,alignValues,varargin)
    defaults = cell(0,3);
    defaults(end+1,:) = {'rescaleShifts','float',1};
    pars = ParseVariableArguments(varargin,defaults,mfilename);
    if isfield(alignValues,'rescale')

        sc =alignValues.rescale; % 1-(1-alignValues.rescale)*pars.rescaleShifts;
        if pars.invert
           sc = 1/sc; 
        end
        S = [sc   0   0;
             0   sc  0;
             0   0   1];
        im = imwarp(im,affine2d(S) ,'OutputView',imref2d(size(im))  );
    end
end

function im = ImTrans(im,alignValues,varargin)
    defaults = cell(0,3);
    defaults(end+1,:) = {'rescaleShifts','float',1};
    pars = ParseVariableArguments(varargin,defaults,mfilename);
    if pars.invert
       pars.rescaleShifts = -pars.rescaleShifts;
    end
    if isfield(alignValues,'xshift') && isfield(alignValues,'yshift')
        im = ImageTranslate(im,[pars.rescaleShifts*alignValues.xshift,...
                                   pars.rescaleShifts*alignValues.yshift]);
%         im = TranslateImage(im,round(pars.rescaleShifts*alignValues.xshift),...
%             round(pars.rescaleShifts*alignValues.yshift));
%  imageTranslate handles 
    end
end

function im = ImRot(im,alignValues,varargin)
    defaults = cell(0,3);
    defaults(end+1,:) = {'rescaleShifts','float',1};
    pars = ParseVariableArguments(varargin,defaults,mfilename);
    if pars.invert
       alignValues.theta = -alignValues.theta; 
    end
    if isfield(alignValues,'theta')
        im = imrotate(im,alignValues.theta,'bilinear','crop'); 
    end
end