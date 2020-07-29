function writeFile = OverwriteDlgBox(filename)
% if the file exists, open a dialog box to confirm if it should be
% overwritten or not.

writeFile = true;
if exist(filename) %#ok<EXIST>
   choice = questdlg(...
       ['File ',filename,' exists.'],... % prompt in the box 
       'Warning',...                              % tile of the box
       'Overwrite','Cancel',...                   % button names
       'Cancel');                                 % default button
    switch choice
        case 'Overwrite'
            writeFile = true;
        case 'Cancel'
            writeFile = false;
    end
end