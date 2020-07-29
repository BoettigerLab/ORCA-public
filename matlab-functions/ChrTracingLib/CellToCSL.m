function listOut = CellToCSL(cellIn)
% Convert cell-array of stings to comma-seperated list
skip = false;

[h,w] = size(cellIn);
if h == 1
    cellIn = cellIn';
    h = w;
elseif isempty(cellIn)
    listOut = '';
    skip = true;
elseif w ~= 1 && h~=1
    error('cellIn must be 1D-cell array');
end

if ~skip
    temp = cat(2,cellIn,repmat({','},h,1));
    temp = reshape(temp',1,2*h);
    listOut = cat(2,temp{1:end-1});
end
% cellIn = strsplit(listOut,',')
    