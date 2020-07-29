function structOut = AddStructures(struct1,struct2,varargin)
% Takes twp structures which contain the same N fields all of a common 
% length.  
% Returns a single structure which has N fields. Each field contains all of
% the entries of structure 1, then all the entries of structure 2
% 
%% Example
% mlist1.x = rand(10,1); mlist1.y = rand(10,1);
% mlist2.x = rand(10,1); mlist2.y = rand(10,1);
% mlist = AddStructures(mlist1,mlist2)
% mlist.x is 20x1 and mlist.y is 20x1 combining the values of mlist 1 and
% mlist 2. 

% Added pad.  If second array is smaller, pad missing data with NaN

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean',true};
defaults(end+1,:) = {'catdim', 'positive',1};
defaults(end+1,:) = {'pad', 'boolean',false};  % needs testing

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'required: struct1,struct2');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


    
%-------------------------------------------------------------------------
%% Main function
%-------------------------------------------------------------------------   


for f=fieldnames(struct1)'
    try
        structOut.(f{1})=cat(parameters.catdim,struct1.(f{1}),struct2.(f{1}));
    catch
        if parameters.verbose
           disp(['field ',f{1},' could not be concatinated']);  
        end
        if parameters.pad
            try
                [h1,w1,d1] = size(struct1.(f{1}));
                [h2,w2,d2] = size(struct2.(f{1}));
                s2 = padarray( struct2.(f{1}), [h1-h2,w1-w2,d1-d2],nan,'post');
                structOut.(f{1})=cat(parameters.catdim,struct1.(f{1}),s2);
            catch
                structOut.(f{1})={struct1.(f{1}),struct2.(f{1})};    
            end
        else
            structOut.(f{1})={struct1.(f{1}),struct2.(f{1})};
        end
    end
end

