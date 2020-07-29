%% h = HEXSCATTER( x, y, ... )
% Gordon Bean, February 2014
% A scatter-plot substitute - generate a density plot using hexagonal
% patches.
%
% Syntax
% hexscatter(xdata, ydata)
% hexscatter(xdata, ydata, 'Name', Value, ...)
% h = hexscatter(...)
%
% Description
% hexscatter(xdata, ydata) creates a density plot of the ydata versus the
% xdata using hexagonal tiles. xdata and ydata should be vectors. NaN
% values (and their corresponding values in the other vector) are ignored.
%
% hexscatter(xdata, ydata, 'Name', Value, ...) accepts name-value pairs of
% arguments from the following list (defaults in {}):
%  'xlim' { [min(xdadta(:) max(xdata(:))] } - a 2-element vector containing
%  the lower and upper bounds of the 2nd dimension of the grid.
%  'ylim' { [min(ydadta(:) max(ydata(:))] } - a 2-element vector containing
%  the lower and upper bounds of the 1st dimension of the grid.
%  'res' { 50 } - the resolution, or number of bins in each dimension. The
%  total number of bins will be the resolution squared.
%  'drawEdges' { false } - if true, edges are drawn around each hexagonal
%  patch.
%  'showZeros' { false } - if true, bins with 0 counts are shaded; if
%  false, only bins with non-zero counts are colored. 
% 
% h = hexscatter( ... ) returns the object handle to the patch object
% created.
% 
% Examples
% hexscatter(rand(2000,1), rand(2000,1))
%
% hexscatter(rand(2000,1), rand(2000,1), 'res', 90)
%
% Also available in the Bean Matlab Toolkit:
% https://github.com/brazilbean/bean-matlab-toolkit
%
%
% Copyright (c) 2014, Gordon
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function h = hexscatter( xdata, ydata, varargin )
    xdata(isinf(xdata)) = nan;
    ydata(isinf(ydata)) = nan;
    params = default_param( varargin, ...
        'xlim', [min(xdata(:)) max(xdata(:))], ...
        'ylim', [min(ydata(:)) max(ydata(:))], ...
        'res', 50, ...
        'drawEdges', false, ...
        'showZeros', false);
    
    if params.drawedges
        ec = 'flat';
    else
        ec = 'none';
    end
    
    %% Determine grid
    xl = params.xlim;
    yl = params.ylim;
    
    xbins = linspace(xl(1), xl(2), params.res);
    ybins = linspace(yl(1), yl(2), params.res);
    dy = diff(ybins([1 2]))*0.5;
    
    [X, Y] = meshgrid(xbins, ybins);
    n = size(X,1);
    Y(:,1:fix(end/2)*2) = ...
        Y(:,1:fix(end/2)*2) + repmat([0 dy],[n,fix(n/2)]);

    %% Map points to boxes
    nix = isnan(xdata) | isnan(ydata);
    xdata = xdata(~nix);
    ydata = ydata(~nix);
    
    % Which pair of columns?
    dx = diff(xbins([1 2]));
    foox = floor((xdata - xbins(1)) ./ dx)+1;
    foox(foox > length(xbins)) = length(xbins);
    
    % Which pair of rows?
    % Use the first row, which starts without an offset, as the standard
    fooy = floor((ydata - ybins(1)) ./ diff(ybins([1 2])))+1;
    fooy(fooy > length(ybins)) = length(ybins);
    
    % Which orientation
    orientation = mod(foox,2) == 1;

    % Map points to boxes
    foo = [xdata - xbins(foox)', ydata - ybins(fooy)'];

    % Which layer
    layer = foo(:,2) > dy;

    % Convert to block B format
    toflip = layer == orientation;
    foo(toflip,1) = dx - foo(toflip,1);

    foo(layer==1,2) = foo(layer==1,2) - dy;

    % Find closest corner
    dist = sqrt(sum(foo.^2,2));
    dist2 = sqrt(sum(bsxfun(@minus, [dx dy], foo).^2, 2));

    topright = dist > dist2;

    %% Map corners back to bins
    % Which x bin?
    x = foox + ~(orientation == (layer == topright));
    x(x > length(xbins)) = length(xbins);
    
    % Which y bin?
    y = fooy + (layer & topright);
    y(y > length(ybins)) = length(ybins);
    
    ii = sub2ind(size(X), y, x);

    %% Determine counts
    counts = sum(bsxfun(@eq, ii, 1:numel(X)),1);

    newplot;
    xscale = diff(xbins([1 2]))*2/3;
    yscale = diff(ybins([1 2]))*2/3;
    theta = 0 : 60 : 360;
    x = bsxfun(@plus, X(:), cosd(theta)*xscale)';
    y = bsxfun(@plus, Y(:), sind(theta)*yscale)';
    
    if params.showzeros
        h = patch(x, y, counts, 'edgeColor', ec);
    else
        jj = counts > 0;
        h = patch(x(:,jj), y(:,jj), counts(jj), 'edgeColor', ec);
    end
    
    if nargout == 0
        clear h;
    end
    
    %% Function: default_param
    % Gordon Bean, March 2012
    % Copied from https://github.com/brazilbean/bean-matlab-toolkit
    function params = default_param( params, varargin )
        if (iscell(params))
            params = get_params(params{:});
        end
        defaults = get_params(varargin{:});

        for f = fieldnames(defaults)'
            field = f{:};
            if (~isfield( params, lower(field) ))
                params.(lower(field)) = defaults.(field);
            end
        end
    end

    %% Function: get_params - return a struct of key-value pairs
    % Gordon Bean, January 2012
    %
    % Usage
    % params = get_params( ... )
    %
    % Used to parse key-value pairs in varargin - returns a struct.
    % Converts all keys to lower case.
    %
    % Copied from https://github.com/brazilbean/bean-matlab-toolkit
    function params = get_params( varargin )
        params = struct;

        nn = length(varargin);
        if (mod(nn,2) ~= 0)
            error('Uneven number of parameters and values in list.');
        end

        tmp = reshape(varargin, [2 nn/2]);
        for kk = 1 : size(tmp,2)
            params.(lower(tmp{1,kk})) = tmp{2,kk};
        end
    end
end