% read nnrd
% fldr_nrrd = '/home/jingjing/Documents/code2022/IR_vessel/fromDjazia/new_merged';
% fn_nrrd = [fldr_nrrd '/merged.nrrd'];

fldr_nrrd = '/home/jingjing/Documents/code2022/IR_vessel/fromDjazia/last_chance'
fn_nrrd = [fldr_nrrd '/this_one.nrrd'];
addpath('/home/jingjing/Documents/SoftwarePackages/nrrdread')

[X, META]  = nrrdread(fn_nrrd);
%%
figure,
for ii = 1:136
    imagesc(squeeze(X(:,:,ii)))
    title([num2str(ii)])
    colorbar
    caxis([0 5])
    pause(0.3)
end

%% visualize single slice
ii =109 
figure
imagesc(squeeze(X(:,:,ii)))
title([num2str(ii)])
    colorbar
    caxis([0 5])
%% save to tiff
[n1 n2 n3] = size(X);
V = uint16(X);
fldr_tif = './tiff_0517/'
mkdir(fldr_tif )
name_tif = 'img_';
 
n_end = 120;
for ii = 1:n_end
I = squeeze(V(:,:,ii));
% title(num2str(ii))
% pause(0.3)
% imwrite(I, ...
%     [fldr_tif '/' name_tif num2str(n_end-ii+1) '.tif'])
imwrite(I, ...
    [fldr_tif '/' name_tif num2str(ii) '.tif'])
end
