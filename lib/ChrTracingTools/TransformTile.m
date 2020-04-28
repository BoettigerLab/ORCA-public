function tilePlaced = TransformTile(refImage,imageTile,ul,alignValues,varargin)
% tilePlaced = TransformTile(refImage,imageTile,ul,tform)
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[h_m,w_m,z_m] = size(refImage);
[h_i,w_i,~] = size(imageTile);
tilePlaced = zeros(h_m,w_m,z_m); % ,class(refImage)

% if initial tile sits outside of refImage, just trim the tile.
  % it is possible this is problematic for image rescaling, that needs
  % thought. 
if isempty(alignValues)
   alignValues.rescale = 1; 
end
  
  
try
    rescaleShiftULx = w_i/2 - alignValues.rescale*w_i/2;
    rescaleShiftULy = h_i/2 - alignValues.rescale*h_i/2;
    ul = ul - [rescaleShiftULx,rescaleShiftULy];
    ul = round(ul);

    if ul(1) <= 0
       if pars.verbose
           disp(['imageTile initial position is outside left of refImage. Trimming image tile by ', num2str(ul(1))]);
       end
       imageTile = imageTile(:,-ul(1)+1:end,:);
       ul(1) = 1;
    end
    if ul(2) <= 0
       if pars.verbose
           disp(['imageTile initial position is outside top of refImage. Trimming image tile by ',num2str(ul(2))]);
       end
        imageTile = imageTile(-ul(2)+1:end,:,:);
        ul(2) = 1;
    end

    [h_i,w_i,~] = size(imageTile);
    if ~isempty(alignValues)
        imageTile = ScaleRotateShift(imageTile,alignValues);
    end
    tilePlaced(ul(2):ul(2)+h_i-1,ul(1):ul(1)+w_i-1,:) = imageTile;
    % tilePlaced = ScaleRotateShift(tilePlaced,alignValues);
    % swapped order
    [h_2,w_2,~] = size(tilePlaced);
    if h_2~=h_m || w_2 ~= w_m
        tilePlaced = tilePlaced(1:h_m,1:w_m,:); % enforce sizing 
        if pars.verbose
            disp(['imageTile final position is outside right or bottom of refImage. Trimming tile by ',num2str([h_2,w_2]-[h_m,w_m])]);
        end
    end
catch er
    disp(er.getReport)
    stopNow = true;
    disp(stopNow);
end
