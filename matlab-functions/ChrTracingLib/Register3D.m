function [outputImage,shifts] = Register3D(refImage,inputImage,varargin)
% [outputImage,shifts] = Register3D(refImage,inputImage)
% 
% i3 = Register3D(i1,i2,'upsample',4,'showplots',false)
% 
% Alistair Boettiger
% Aug 4, 2017
% CC BY NC

defaults = cell(0,3);
defaults(end+1,:) = {'upsample','positive',4};
defaults(end+1,:) = {'upsampleZ','positive',4};
defaults(end+1,:) = {'showplots','boolean',false};
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'figShowAlign','freeType',12};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'maxShiftXY','positive',inf};
defaults(end+1,:) = {'maxShiftZ','positive',inf};
defaults(end+1,:) = {'center','nonnegative',0};
defaults(end+1,:) = {'threshold','nonnegative',.6};
defaults(end+1,:) = {'use3D','boolean',false};
defaults(end+1,:) = {'returnImage','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% allow the feature to be turned off completely
if pars.maxShiftXY == 0 && pars.maxShiftZ == 0
    outputImage = inputImage;
else

    % start main function
    [nRows,nCols,nStks] = size(refImage);
    [iRows,iCols,iStks] = size(inputImage);
    if pars.center==0
        pars.center = [nRows/2,nCols/2,nStks/2]+.5;
    end

    % protect against over trimming the box
    bXY = max(4,pars.maxShiftXY);
    bZ = max(4,pars.maxShiftZ);
    
    xi = max(round(pars.center(1)-bXY),1);
    xe = min([round(pars.center(1)+bXY),nCols,iCols]);
    yi = max(round(pars.center(2)-bXY),1);
    ye = min([round(pars.center(2)+bXY),nRows,iRows]);
    zi = max(round(pars.center(3)-bZ), 1);
    ze = min([round(pars.center(3)+bZ), nStks, iStks]);

    i1 = refImage(yi:ye,xi:xe,zi:ze);
    i2 = inputImage(yi:ye,xi:xe,zi:ze);
    
    [rowsi,colsi,zsi] = size(i1);
    rowsI = round(rowsi*pars.upsample);
    colsI = round(colsi*pars.upsample);
    zsI = round(zsi*pars.upsampleZ);
    im1 = imresize3(i1,[rowsI,colsI,zsI]);
    % im1= imresize3(i1,pars.upsample); % try imresize3, change all 3 dimensions 

    bx =  [i1(1,:,:), i1(end,:,:)];
    by =  [i1(:,1,:), i1(:,end,:)];
    bz =  [i1(:,:,1), i1(:,:,end)];
    edge1 = quantile( [bx(:) ; by(:); bz(:)],.9); % .8

    bx =  [i2(1,:,:), i2(end,:,:)];
    by =  [i2(:,1,:), i2(:,end,:)];
    bz =  [i2(:,:,1), i2(:,:,end)];
    edge2 = quantile( [bx(:) ; by(:); bz(:)],.9); % .8

    % Get a small volume as template
    [rowsi,colsi,zsi] = size(i2);
    rowsI = round(rowsi*pars.upsample);
    colsI = round(colsi*pars.upsample);
    zsI = round(zsi*pars.upsampleZ);
    im2 = imresize3(i2,[rowsI,colsI,zsI]);
    
    % im2= imresize3(i2,pars.upsample);  % try imresize3, change all 3 dimensions 
    %im2 = ImageTranslate(im2,[0,0]);
    

    im1o = im1; % save original for plotting
    im2o = im2; % save original for plotting

    
    % KEY: send edge values to zero, because we will pad with zeros when
    % doing the correlation based alignment.  
    if pars.threshold ~= 0
        im1 = im1 - 1*edge1;
        im2 = im2 - 1*edge2;
    else
        theta = quantile(im1(:), pars.threshold);
        im1(im1 < theta) = 0;
        theta = quantile(im2(:), pars.threshold);
        im2(im2<theta) = 0;
    end

    %---------------------------% 
    % Full 3D
    %----------------------------
    if pars.use3D
        % Calculate SDD between template and image
        I_SSD=template_matching(im2,im1);
        % Find maximum correspondence
        [shifts.score,indmax] = max(I_SSD(:));
        [y,x,z]=ind2sub(size(I_SSD),indmax);
        % [y,x,z]=ind2sub(size(I_SSD),find(I_SSD==max(I_SSD(:))));
        [h,w,d] = size(im1);
        xshift = x-w/2;
        yshift = y-h/2;
        zshift = z-d/2;
        %--------------------------%    
        shifts.xshift = xshift/pars.upsample;
        shifts.yshift = yshift/pars.upsample;
        shifts.zshift = zshift/pars.upsampleZ; % updated with imresize3   
        maxedOut = shifts.xshift >= pars.maxShiftXY ||...
                   shifts.yshift >= pars.maxShiftXY ||...
                   shifts.zshift >= pars.maxShiftZ;
    else
        maxedOut = true;
        shifts.score = 0;
    end
   
   if maxedOut || shifts.score < .99    
       % Try 2D first:
        im1_xy = max(im1,[],3);
        im2_xy = max(im2,[],3);
        if pars.showExtraPlots; figure(13); clf; end
        alignData_xy = CorrAlignFast(im1_xy,im2_xy,'showplot',pars.showExtraPlots,'gradMax',false);
        im1_yz = max(permute(im1,[3,2,1]),[],3);
        im2_yz = max(permute(im2,[3,2,1]),[],3);
        if pars.showExtraPlots; figure(14); clf; end
        alignData_yz = CorrAlignFast(im1_yz,im2_yz,'showplot',pars.showExtraPlots,'gradMax',false);
        xshift = mean([alignData_xy.xshift,alignData_yz.xshift]);
        yshift = alignData_xy.yshift; %  mean([,alignData_yz.xshift]);
        zshift = alignData_yz.yshift;
        % rescale
        shifts.xshift = xshift/pars.upsample;
        shifts.yshift = yshift/pars.upsample;
        shifts.zshift = zshift/pars.upsampleZ; % updated with imresize3
        shifts.score = 1;      
       % test if maxed out
         maxedOut = shifts.xshift >= pars.maxShiftXY ||...
               shifts.yshift >= pars.maxShiftXY ||...
               shifts.zshift >= pars.maxShiftZ;
       % don't accept larger than max shifts
       if maxedOut      
          if pars.verbose
             disp('registration failed, requested shift >= max shift.');
          end
          xshift = 0; yshift = 0; zshift = 0;
          shifts.xshift = 0; shifts.yshift =0; shifts.zshift = 0;
          shifts.score = 0;
       end         
   end
           
    if pars.verbose
        disp(['computed ',...
            ' xshift ', num2str(shifts.xshift),...
            ' yshift ',num2str(shifts.yshift),...
            ' zshift ',num2str(shifts.zshift),...
            ' score: ',num2str(shifts.score)]);
    end
    
    if pars.returnImage || pars.showplots
        [rowsi,colsi,zsi] = size(inputImage);
        rowsI = round(rowsi*pars.upsample);
        colsI = round(colsi*pars.upsample);
        zsI = round(zsi*pars.upsampleZ);
        imTemp = imresize3(inputImage,[rowsI,colsI,zsI]);
        % imTemp = imresize3(inputImage,pars.upsample);
        temp = TranslateImage(imTemp,...
         xshift,yshift,'zshift',zshift);
     
        [rowsi,colsi,zsi] = size(temp);
        rowsI = round(rowsi/pars.upsample);
        colsI = round(colsi/pars.upsample);
        zsI = round(zsi/pars.upsampleZ);
        i3 = imresize3(temp,[rowsI,colsI,zsI]);
     
        % i3 = imresize3(temp,1/pars.upsample); 
        i3(i3==0) = edge2;
        outputImage = i3;
    else
       outputImage = []; 
    end

    if pars.showplots
        if pars.showExtraPlots
            figure(11); clf;
            subplot(2,2,1); Ncolor(IncreaseContrast( cat(3, max(refImage,[],3), max(inputImage,[],3)) ) ); xlabel('x'); ylabel('y'); title('ref vs input xy')
            hold on; plot(pars.center(1),pars.center(2),'b+');
            subplot(2,2,2); Ncolor(IncreaseContrast( cat(3, max(permute(refImage,[3,2,1]),[],3), max(permute(inputImage,[3,2,1]),[],3) ) )); xlabel('y'); ylabel('z'); title('ref vs input yz');
            hold on; plot(pars.center(2),pars.center(3),'b+');
            subplot(2,2,3); Ncolor(IncreaseContrast( cat(3, max(refImage,[],3), max(outputImage,[],3)) ));  title('ref vs output xy');
            hold on; plot(pars.center(1),pars.center(2),'b+');
            subplot(2,2,4); Ncolor(IncreaseContrast( cat(3, max(permute(refImage,[3,2,1]),[],3), max(permute(outputImage,[3,2,1]),[],3) ) )); xlabel('y'); ylabel('z'); title('ref vs output yz')
            hold on; plot(pars.center(2),pars.center(3),'b+');
        end
        if pars.figShowAlign
            figure(pars.figShowAlign); clf;
            subplot(2,2,1); Ncolor(IncreaseContrast( cat(3, max(im1o,[],3), max(im2o,[],3)) ) );  title('upsample ref vs output xy')
            subplot(2,2,2); Ncolor(IncreaseContrast( cat(3, max(permute(im1o,[3,2,1]),[],3), max(permute(im2o,[3,2,1]),[],3) ) )); xlabel('y'); ylabel('z');
            im3 = TranslateImage(im2,xshift,yshift ,'zshift',zshift);
            subplot(2,2,3); Ncolor(IncreaseContrast( cat(3, max(im1,[],3), max(im3,[],3)) ));
            subplot(2,2,4); Ncolor(IncreaseContrast( cat(3, max(permute(im1,[3,2,1]),[],3), max(permute(im3,[3,2,1]),[],3) ) )); xlabel('y'); ylabel('z');
        end
    end
end

% 