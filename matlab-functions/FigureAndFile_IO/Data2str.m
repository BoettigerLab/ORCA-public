function s = Data2str(x,precision)
%INT2STR Convert integer to string.
%   S = INT2STR(X) rounds the elements of the matrix X to integers and 
%   converts the result into a string matrix.
%   Return NaN and Inf elements as strings 'NaN' and 'Inf', respectively.
%
%   See also NUM2STR, SPRINTF, FPRINTF, MAT2STR.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.20.4.11 $  $Date: 2010/11/22 02:46:49 $

% only work with real portion of x
x = real(x);

%precision = 3;

% reduce to indicated precision
digits = round(log10(x));
x = round(x/10^(digits - precision+1))*10^(digits - precision+1);


% create a copy of x to use to calculate maximum width in digits
widthCopy = x;
if isfloat(x)
    x = 0+round(x); %remove negative zero
    % replace Inf and NaN with a number of equivalent length for width
    % calcultion
    widthCopy(~isfinite(widthCopy)) = 314;
    formatConversion = '.0f';
elseif isa(x, 'uint64')
    formatConversion = 'lu';
else
    formatConversion = 'ld';
end

if isempty(x)
    s = '';
elseif isscalar(x)
    s = sprintf(['%', formatConversion], x);
else
    % determine the variable text field width quantity
    widthMax = double(max(abs(widthCopy(:))));
    if widthMax == 0
        width = 3;
    else
        width = floor(log10(widthMax)) + 3;
    end

    format = sprintf('%%%d%s', width, formatConversion);

    [rows, cols] = size(x);
    s = char(zeros(rows, width*cols));
    for row = 1:rows
        % use vectorized version of sprintf for each row
        s(row,:) = sprintf(format, x(row,:));
    end

    % trim leading spaces from string array within constraints of rectangularity.
    s = strtrim(s);
end

