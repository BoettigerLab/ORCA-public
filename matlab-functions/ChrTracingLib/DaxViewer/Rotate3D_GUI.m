function figureHandle = Rotate3D_GUI(im3D,varargin)
% Interactively rotate a 3D image, show a max

defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);
if isempty(pars.figHandle)
    figureHandle = figure(200); clf;
else
    try
        figureHandle = figure(pars.figHandle);
    catch er
       disp(er.getReport);
       error('passed figure handle is not a valid fig handle!') 
    end
end
imagesc(max(im3D,[],3));
userData.image = im3D;
userData.clicked = false;
userData.startPos = [0,0,0];
userData.rot = [0,0]; % initialize
figureHandle.UserData = userData;
figureHandle.WindowButtonMotionFcn = @motionCallback;
figureHandle.WindowButtonDownFcn = @downCallback;
figureHandle.WindowButtonUpFcn = @upCallback;

function downCallback(obj,evd) %#ok<*INUSD>
cp=get(gca,'CurrentPoint');
obj.UserData.clicked = true;
obj.UserData.startPos = cp(1,1:2);

function upCallback(obj,evd)
obj.UserData.clicked = false;
new3D = RotateMatrix3D(obj.UserData.image,obj.UserData.rot(1),obj.UserData.rot(2));
imagesc(max(new3D,[],3));

function motionCallback(obj,evd)
if obj.UserData.clicked
    cp=get(gca,'CurrentPoint');
    mov = obj.UserData.startPos - cp(1,1:2);    
    tx = rem(mov(1)*10,360);
    ty = rem(mov(2)*10,360);
    obj.UserData.rot = [tx,ty];
    new3D = RotateMatrix3D(obj.UserData.image,tx,ty);
    imagesc(max(new3D,[],3));
end
