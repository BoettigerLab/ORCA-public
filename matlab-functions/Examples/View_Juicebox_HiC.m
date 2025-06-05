%% Example load Hi-C data

% the HiC data has been exported from juicebox as a .txt file using the
% export function in juicebox.  It's not necessary to get the window to
% match the domain of interest perfectly, as long as the domain is visible
% (on both axes!) prior to export. Note it should be saved as .txt. 

mapRes = 25E3;
HoxA_domain3Mb = 'chr6:50,590,001-52,780,076'; % mm10  % 
jboxFolder = 'U:\GenomeData\JuiceboxExport\';

[hic_hoxA_cc,coords] = ReadJuiceboxMatrix([jboxFolder,'mESC_HoxA_CTCFcntrl_25kb.txt'],...   
                    'locus',HoxA_domain3Mb,'mapRes',mapRes,'displayRes',mapRes); 

[hic_hoxA_ctcf,coords] = ReadJuiceboxMatrix([jboxFolder,'mESC_HoxA_CTCFdeg_25kb.txt'],...  
                    'locus',HoxA_domain3Mb,'mapRes',mapRes,'displayRes',mapRes);

% you may chose to downsample the resolution (e.g. displayRes might be
% 2*mapRes), but it's too late to upsample.  

% display data
figure(1); clf; 
subplot(1,2,1); imagesc((hic_hoxA_cc)); 
clim([0,.2*max(hic_hoxA_cc(:))]);  colorbar; title('untreated');
subplot(1,2,2); imagesc((hic_hoxA_ctcf));  
clim([0,.2*max(hic_hoxA_ctcf(:))]);  colorbar;   title('CTCF degraded')
GetColorMap('whiteToRed');

% subtract raw maps
figure(2); clf; imagesc(hic_hoxA_cc-hic_hoxA_ctcf);
GetColorMap('RedWhiteBlueSat'); colorbar; clim([-4000,4000]);
title('difference (UT - Aux)')

% balance maps and subtract
norm_cc = hic_hoxA_cc/mean(hic_hoxA_cc(:));
norm_ctcf = hic_hoxA_ctcf/mean(hic_hoxA_ctcf(:));
figure(2); clf; imagesc(norm_cc-norm_ctcf);
GetColorMap('RedWhiteBlue'); colorbar; clim([-3,3]);
title('balanced UT - Aux')

% balance maps and log ratio
figure(3); clf; imagesc(log2(norm_cc./norm_ctcf));
GetColorMap('RedWhiteBlue'); colorbar; clim([-3,3]);
title('balanced log2 ratio (UT/Aux)')