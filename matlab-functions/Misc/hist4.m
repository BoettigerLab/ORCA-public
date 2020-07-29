
function M = hist4(x,y,z,varargin)
% M = hist4(x,y,z)
% create a 3D density plot (3d histogram) from 3D data (x,y,z).  
% bins 
%-------------------------------------------------------------------------
% Required Inputs:
% x,y,z -- vector coordinates of the 3D data.  
%
%-------------------------------------------------------------------------
% Outputs:
%  M is a HxWxN matrix, where each voxel contains the number of points in x
%  y,z which fall into it.  
%
%-------------------------------------------------------------------------
% Optional inputs: ('name' / datatype / default)
% 'bins' / scalar or vector / [100,100,100]
%                      -- number of bins in each dimension.  If a scalar,
%                      all dimensions will have the given number of bins.
% 'datarange' / cell / min to max of data
%                       -- cell of two dimensional vectors [xi-min,xi-max].
%                       Changes the min and max range to include in the
%                       3D histogram for any dimension. 
%
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger.alistair@gmail.com
% February 19th, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License 3.0 CC BY  
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
% Default Parameters
%-------------------------------------------------------------------------
bins = [100,100,100];
datatype = 'uint16';
datarange = cell(1,3); 

%-------------------------------------------------------------------------
% Parse variable input
%-------------------------------------------------------------------------

if nargin > 3
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'bins'
                bins = CheckParameter(parameterValue, 'positive', 'bins');
            case 'datatype'
                datatype = CheckParameter(parameterValue, 'string', 'datatype');
            case 'datarange'
                datarange = CheckParameter(parameterValue, 'cell', 'datarange');
            otherwise
                error(['The parameter ''', parameterName,...
                    ''' is not recognized by the function, ''',...
                    mfilename '''.' '  See help ' mfilename]);
        end
    end
end


%-------------------------------------------------------------------------
%% Main Function
%-------------------------------------------------------------------------   

if length(bins) == 1
    bins = repmat(bins,1,3);
end

[rx,ry,rz] = datarange{:};

xbins = bins(1);
ybins = bins(2);
zbins = bins(3);

if isempty(rz)
    rz = [min(z),max(z)]; % range of z
end
if isempty(ry)
    ry = [min(y),max(y)];
end
if isempty(rx)
    rx = [min(x),max(x)];
end

Zs = linspace(rz(1),rz(2),zbins);
Ys = linspace(ry(1),ry(2),ybins);
Xs = linspace(rx(1),rx(2),xbins);
Zs = [Zs,inf];

M = zeros(ybins,xbins,zbins,datatype); 
for i=1:zbins % i=6
    inplane = z>Zs(i) & z<Zs(i+1);
    yi = y(inplane);
    xi = x(inplane);
    if ~isempty(yi) & ~isempty(xi)
        M(:,:,i) = hist3([yi,xi],{Ys,Xs});
    end
end


%    [~,~,zs] = size(M);
%         figure(10);
%         k=0;
%         for j=1:zs
%             k=k+1;
%             subplot(6,6,k); imagesc(M(:,:,j));
%         end
%         

