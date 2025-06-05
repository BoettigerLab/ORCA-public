function r = CorCoef(x,y,varargin)
% because it's annoying to always need two lines to do this
%   we just want a single correlation coefficient, not a 2x2 symmetric
%   matrix with 1s on the diagonal.
% 

try
    if isempty(x) | isempty(y) | length(x) < 2 | length(y) < 2
        r=0;
    else
        if isempty(varargin)
            m = corrcoef(x,y);
        else
            m = corrcoef(x,y,varargin); % needs testing
        end
        r = m(1,2); 
    end
catch er
    warning(er.getReport)
    disp('debug here')
end
