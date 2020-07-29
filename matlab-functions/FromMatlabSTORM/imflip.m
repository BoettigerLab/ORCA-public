function Iout = imflip(Iin,dim)
%  Iout = imflip(Iin,dim)
% 
% flip image along the specified dimension, dim:
%       dim=1 for flip left-right
%       dim=2 for flip up-down
% image may be 2D or 3D
% 
% for multiple flips, stack imflip commands in the desired order
%   imOut = imflip(imflip(imIn,1),2)

% updated 04/28/11 to use single channel as well as multichannel data.


chns = size(Iin,3);

if chns == 1
    if dim == 1
        Iout = fliplr(Iin);
    elseif dim == 2
        Iout = flipud(Iin);
    elseif dim == 0
        Iout = Iin; 
    end
else

    intype = class(Iin);
Iout = zeros(size(Iin),intype);
if dim == 1
    I1 = fliplr(Iin(:,:,1));
    I2 = fliplr(Iin(:,:,2));
    I3 = fliplr(Iin(:,:,3));
    Iout(:,:,1) = I1;
    Iout(:,:,2) = I2;
    Iout(:,:,3) = I3;
elseif dim == 2
    I1 = flipud(Iin(:,:,1));
    I2 = flipud(Iin(:,:,2));
    I3 = flipud(Iin(:,:,3));
    Iout(:,:,1) = I1;
    Iout(:,:,2) = I2;
    Iout(:,:,3) = I3;
elseif dim == 0
    Iout = Iin; 
else
    disp('error: expected 2nd entry as dimension = 1 or 2 or 0 for no flip'); 
end

end