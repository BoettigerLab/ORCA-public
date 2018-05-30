function struct3 = JoinStructures(struct1,struct2,varargin)
% take two structures, return a third that has all the entries and field
% names of the first two.  
% 
%% Example
% clear struct1 struct2;
% struct1.a = 1;
% struct1.z = 2;
% struct1.A = 'Hello';
% struct2.b = 3;
% struct2.z = 6;
% struct2.B = 'Goodbye';
% struct3 = JoinStructures(struct1,struct2,'conflict','keepSecond');
% struct3 = 
%   struct with fields: 
%     a: 1
%     A: 'Hello'
%     b: 3
%     z: 6
%     B: 'Goodbye'

defaults = cell(0,3);
defaults(end+1,:) = {'conflict',{'keepFirst','keepSecond'},'keepSecond'};
pars = ParseVariableArguments(varargin,defaults,mfilename); 

% 
names1 = fieldnames(struct1); 
names2 = fieldnames(struct2); 

if strcmp(pars.conflict,'keepSecond')
    [v,i] = intersect(names1,names2);
    names1(i) = [];
    for n=1:length(v)
        struct1 = rmfield(struct1,v{n});
    end
else
    [v,i] = intersect(names2,names1);
    names2(i) = [];
    for n=1:length(i)
        struct2 = rmfield(struct2,v{n});
    end
end

names = [names1; names2];
struct3 = cell2struct([struct2cell(struct1); struct2cell(struct2)], names, 1);

