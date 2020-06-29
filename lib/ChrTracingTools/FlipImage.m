function imageOut = FlipImage(imageIn,flipdim)
% 1 = flip horizontal
% 2 = flip vertical
% 3 = transpose

chns = size(imageIn,3);
imageOut = imageIn;

if flipdim == 1 
    for i=1:chns
           imageOut(:,:,i) = fliplr(imageIn(:,:,i));
    end
elseif flipdim == 2
   for i=1:chns
       imageOut(:,:,i) = flipud(imageIn(:,:,i));
   end
elseif flipdim == 3
   for i=1:chns
       imageOut(:,:,i) = imageIn(:,:,i)';
   end
end
    