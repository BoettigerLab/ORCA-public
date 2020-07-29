1.  Download following libraries
 This code requires following toolbox/submission available on matlab central
  1.1. Toolbox Wavelets by Gabriel Peyre available at http://www.mathworks.com/matlabcentral/fileexchange/5104-toolbox-wavelets'
  1.2. k-means clustering available at http://in.mathworks.com/matlabcentral/fileexchange/24616-kmeans-clustering/content/kmeans/kmeans.m . Rename kmeans.m as litekmeans.m to distinguish it from matlab inbuilt kmeans.

2. Run the program WaveletBasedImageSegmentation.m

3. To test on your image, change the name of image file in "I=imread('syn_new_9.bmp');" at line 22

4. Adjust the Parameters and run the program
      IsGelImage=1;      %tell if it is 2D gel images/ or images with blob like objects for extra processing of binary image with prior knowledge of radius of smallest spot
      IsObjectLighter=0; %set 1 if the objects are lighter than background, otherwise set 0


      if (IsGelImage == 0)
	
	KernelBandwidth_h =1;     
	ContrastThreshold_Th = 3;   % low value means detect low contrast objects, high value means detect high contrast objects
	MorphSmoothRad=1;   
	smallRadiusofSpot=0;  % required in case of 2D gel image.
    
     else
      
	KernelBandwidth_h =6;     
	ContrastThreshold_Th = 6;    
	MorphSmoothRad=3;   %morphological disk radius set to smaller than small radius of spot.
	smallRadiusofSpot=5;  % required in case of 2D gel image.

      end
5. For high contrast object in noise free image 
	KernelBandwidth_h =1;     
	ContrastThreshold_Th = 45;   % high value means detect high contrast objects try different values then fix it
6. For image consisting small size objects, set KernelBandwidth_h =1; 
7. For noisy images/images with inhomogenous background, set KernelBandwidth_h > 1; 