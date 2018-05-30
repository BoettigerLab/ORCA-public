function clrmap = GetColorMap(clrmapName,varargin)
% GetColorMap('hsv')
% GetColorMap('redToWhite');
% GetColorMap('redToWhite',10);
%

global matlabFunctionsPath


% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};


if nargin > 1
    pts = varargin{1};
else
    pts = 256;
end

varin = varargin(2:end);
parameters = ParseVariableArguments(varin,defaults,mfilename);


if ~ischar(clrmapName)
    clrmap = clrmapName;
    return
elseif strcmp(clrmapName,'default')
    try
        clrmap = parula(pts);
    catch
        clrmap = jet(pts);
    end
    return
end

try 
    if strcmp(clrmapName,hsv)
        clrmap = eval([clrmapName, '(',num2str(pts+1),')']);
    else
        clrmap = eval([clrmapName, '(',num2str(pts),')']);
    end
catch

    % Black to white colormaps via the indicated color name;
  switch clrmapName
         
      
        case 'yellow'
        clrmap = hot(pts);
        clrmap = [clrmap(:,1),clrmap(:,1),clrmap(:,2)];
        clrmap(clrmap<0) = 0;

        case 'red'
        clrmap = hot(pts);
        clrmap = [clrmap(:,1),clrmap(:,2),clrmap(:,3)];
        clrmap(clrmap<0) = 0;
        
        case 'blue'
        clrmap = hot(pts);
        clrmap = [clrmap(:,3),clrmap(:,2),clrmap(:,1)];
        clrmap(clrmap<0) = 0;

        case 'green'
        clrmap = hot(pts);
        clrmap = [clrmap(:,3),clrmap(:,1),clrmap(:,2)];
        clrmap(clrmap<0) = 0;

        case 'purple'
        clrmap = hot(pts);
        clrmap = [clrmap(:,1),clrmap(:,3),clrmap(:,1)];
        clrmap(clrmap<0) = 0; 
        
        case 'black'
        clrmap = gray(pts);
        
        case 'cyan'
        clrmap = hot(pts);
        clrmap = [clrmap(:,3),clrmap(:,1),clrmap(:,1)];
        clrmap(clrmap<0) = 0;

      case 'whiteOrangeRed'
        nPts = pts;
        whiteToYellow = zeros(nPts,3);
        yellowToRed  = zeros(nPts,3);
        redToBlack  = zeros(nPts,3);
        redToWhite = zeros(nPts,3); 
        redToYellow= zeros(nPts,3); 
        for n=1:nPts
            yellowToRed(n,:) = [1,(nPts-n+1)/nPts,0];
            redToBlack(n,:) =  [(nPts-n+1)/nPts,0,0];
            redToWhite(n,:) = [1,n/nPts,n/nPts];
            whiteToYellow(n,:) = [1, (nPts-n*.3+.3)/nPts, (nPts-n+1)/nPts];
            redToYellow(n,:)= [1,.7*n/nPts,0];
        end
        clrmap = flipud([redToYellow;flipud(whiteToYellow)]);

      case 'blackCyanOrange'
        nPts = round(pts/2);
        blackToCyan = zeros(nPts,3);
        CyanToOrange  = zeros(nPts,3);
        for n=1:nPts
            blackToCyan(n,:) = [0 n/nPts n/nPts];
            CyanToOrange(n,:) = [ n/nPts, (nPts-(n/2))/nPts, (nPts-n)/nPts];
        end
        clrmap = ([blackToCyan; (CyanToOrange)]);
        
      case 'PurpleWhiteYellow'  
          nPts = round(pts/2);
          purpleToWhite = zeros(nPts,3);
          whiteToYellow = zeros(nPts,3);
          grayMap = zeros(nPts,3);
          for n=1:nPts
              purpleToWhite(n,:) = [1,n/nPts,1];
              whiteToYellow(n,:) = [1,1,(nPts-n+1)/nPts];
              grayMap(n,:) = [n/nPts,n/nPts,n/nPts];
          end
          clrmap = [purpleToWhite*(2/3)+grayMap*(1/3); whiteToYellow*(7/8)+flipud(grayMap)*(1/8)];
          clrmap(clrmap>1) = 1;
          
      case 'RedWhiteBlue'
          nPts = round(pts/2);
          redToWhite = zeros(nPts,3);
          whiteToBlue = zeros(nPts,3);
          for n=1:nPts
              redToWhite(n,:) = [1,n/nPts,n/nPts];
              whiteToBlue(n,:) = [(nPts-n+1)/nPts,(nPts-n+1)/nPts,1];
          end
          clrmap = [redToWhite; whiteToBlue];
      
      case 'RedWhiteBlueK'
          nPts = round(pts/3);
          redToWhite = zeros(nPts,3);
          whiteToBlue = zeros(nPts,3);
          for n=1:nPts
              kRedToRed(n,:) = [.5+n/nPts*.5,0,0];
              kBlueToBlue(n,:) = [0,0,.5+n/nPts*.5];
              redToWhite(n,:) = [1,n/nPts,n/nPts];
              whiteToBlue(n,:) = [(nPts-n+1)/nPts,(nPts-n+1)/nPts,1];
          end
          clrmap = [kRedToRed; redToWhite; whiteToBlue; flipud(kBlueToBlue)];
          
      case 'redToWhiteK'
        nPts = round(pts*2/3);
        redToWhite = zeros(nPts,3);
        for n=1:nPts
          redToWhite(n,:) = [1,((n-1)/(nPts-1)),((n-1)/(nPts-1))];
        end
        nPts = round(pts*1/3);
        kRedToRed = zeros(nPts,3);
        for n=1:nPts
          kRedToRed(n,:) = [.5*(n-1)/(nPts-1)+.5,0,0];
        end
        clrmap =[kRedToRed; redToWhite];   
          
      case 'RedWhiteBlueSat'
          nPts = round(pts/2);
          redToWhite = zeros(nPts,3);
          whiteToBlue = zeros(nPts,3);
          for n=1:nPts
              redToWhite(n,:) = [1,n/nPts,n/nPts];
              whiteToBlue(n,:) = [(nPts-n+1)/nPts,(nPts-n+1)/nPts,1];
          end
          clrmap = [.9 .9 .9; redToWhite; whiteToBlue; .2 .2 .8];

       case 'BlueWhiteRed'
          nPts = round(pts/2);
          redToWhite = zeros(nPts,3);
          whiteToBlue = zeros(nPts,3);
          for n=1:nPts
              redToWhite(n,:) = [1,n/nPts,n/nPts];
              whiteToBlue(n,:) = [(nPts-n+1)/nPts,(nPts-n+1)/nPts,1];
          end
          clrmap = flipud([redToWhite; whiteToBlue]);
          
       case 'BlueWhiteRedSat'
          nPts = round(pts/2);
          redToWhite = zeros(nPts,3);
          whiteToBlue = zeros(nPts,3);
          for n=1:nPts
              redToWhite(n,:) = [1,n/nPts,n/nPts];
              whiteToBlue(n,:) = [(nPts-n+1)/nPts,(nPts-n+1)/nPts,1];
          end
          clrmap = flipud([.8 .2 .2; redToWhite; whiteToBlue; .9 .9 .9]);
          
    case 'redToWhite'
        nPts = pts;
        redToWhite = zeros(nPts,3);
        for n=1:nPts
          redToWhite(n,:) = [1,((n-1)/(nPts-1)),((n-1)/(nPts-1))];
        end
        clrmap = redToWhite;   
        
        

        
    case 'blueToWhite'
        nPts = pts;
        blueToWhite = zeros(nPts,3);
        for n=1:nPts
          blueToWhite(n,:) = [((n-1)/(nPts-1)),((n-1)/(nPts-1)),1];
        end
        clrmap = blueToWhite; 
        
    case 'redToWhiteSat'
        nPts = pts;
        redToWhite = zeros(nPts,3);
        for n=1:nPts
          redToWhite(n,:) = [1,((n-1)/(nPts-1)),((n-1)/(nPts-1))];
        end
        redToWhite = cat(1,[.8 .8 .8],redToWhite,[.96 .96 .96]);
        clrmap = redToWhite;   
        
    case 'redToWhite2'
        nPts = pts;
        redToWhite = zeros(nPts,3);
        for n=1:nPts
          redToWhite(n,:) = [1,((n-1)/(nPts-1))^.5,((n-1)/(nPts-1))^.5];
        end
        clrmap = redToWhite; 
        
    case 'whiteToRed'
        nPts = pts;
        whiteToRed = zeros(nPts,3);
        for n=1:nPts
          whiteToRed(n,:) = [1,((nPts-n)/(nPts-1)),((nPts-n)/(nPts-1))];
        end
        clrmap = whiteToRed;  
        
    
    case 'whiteToRedSat'
        nPts = pts;
        whiteToRed = zeros(nPts,3);
        for n=1:nPts
          whiteToRed(n,:) = [1,((nPts-n)/(nPts-1)),((nPts-n)/(nPts-1))];
        end
        whiteToRed = cat(1,[.8 .8 .8],whiteToRed,[.96 .96 .96]);
        clrmap = whiteToRed;

    case 'hsvCut'
          nPts = pts;
          clrmap = hsv(round(1.1*nPts)); 
          clrmap = clrmap(1:pts,:);
        
    case 'hsvB'
          nPts = pts;
          clrmap = hsv(round(1.2*nPts));
          cut = round(.65*nPts);
          clrmap = [  clrmap(cut:end,:); clrmap(1:cut-round(.2*nPts),:) ];    
  
    case 'hsvG'
          nPts = pts;
          clrmap = hsv(round(1.2*nPts));
          cut = round(.45*nPts);
          clrmap = [  clrmap(cut:end,:); clrmap(1:cut-round(.2*nPts),:) ];    
  
    case 'viridis'
          vir = readtable([matlabFunctionsPath,'Misc\','viridis_rgb_256.txt']);
          clrmap = vir{:,1:3}./256;
          
    case 'viridisFlip'
          vir = readtable([matlabFunctionsPath,'Misc\','viridis_rgb_256.txt']);
          clrmap = flipud(vir{:,1:3}./256);
          
    otherwise
        if parameters.verbose
            warning(['colormap ',clrmapName,' not recognized']);
        end
  end
end

if nargout == 0
        colormap(clrmap); 
end