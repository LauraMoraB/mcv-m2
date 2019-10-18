%% Generate masks
clearvars;
dst = imread('images/surf_800x450.png');
src = imread('images/kira_800x450.png');

[mask_src, poly] = getMask(src);
[im_s2, mask_dst] = alignSource(src, mask_src, dst);

mask_src=logical(mask_src);
mask_dst=logical(mask_dst);
src=double(src);
dst=double(dst);

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

translate_masked_src = get_translated_src(src, mask_src, mask_dst) ;

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    %[drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    [drivingGrad_i, drivingGrad_j] = mixGradients(translate_masked_src(:,:,nC), dst(:,:,nC));
    
    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);
     
    driving_on_dst = zeros(size(dst(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_dst(:));
    
    param.driving = driving_on_dst;
    
    dst1(:,:,nC) = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst, param);
    
end
imshow(dst1/256)
%% LENA AND GIRL
dst = double(imread('images/lena.png'));
src = double(imread('images/girl.png')); % flipped girl, because of the eyes

mask_src=logical(imread('images/mask_src_eyes.png'));
mask_dst=logical(imread('images/mask_dst_eyes.png'));

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;

translate_masked_src = get_translated_src(src, mask_src, mask_dst);

for nC = 1: nChannels

    %TO DO: COMPLETE the ??
    %[drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    [drivingGrad_i, drivingGrad_j] = mixGradients(translate_masked_src(:,:,nC), dst(:,:,nC));

    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);

    driving_on_dst = zeros(size(src(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_dst(:));
    
    param.driving = driving_on_dst;
    dst1(:,:,nC) = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst, param);
end

mask_src=logical(imread('images/mask_src_mouth.png'));
mask_dst=logical(imread('images/mask_dst_mouth.png'));

translate_masked_src = get_translated_src(src, mask_src, mask_dst) ;

for nC = 1: nChannels
    %TO DO: COMPLETE the ??
    %[drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    [drivingGrad_i, drivingGrad_j] = mixGradients(translate_masked_src(:,:,nC), dst(:,:,nC));

    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);
    
    driving_on_dst = zeros(size(src(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_dst(:));
    
    param.driving = driving_on_dst;

    dst1(:,:,nC) = sol_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst, param);
end
imshowpair(dst/256, dst1/256,'montage')
%% 
function [drivingGrad_i, drivingGrad_j] = importGradients(src)
    [srcGrad_i, srcGrad_j] = gradient(src);
    drivingGrad_i = srcGrad_i;
    drivingGrad_j = srcGrad_j;
end

function [drivingGrad_i, drivingGrad_j] = mixGradients(src, dst)
    [srcGrad_i, srcGrad_j] = gradient(src);
    [dstGrad_i, dstGrad_j] = gradient(dst);
    cond = magnitude(srcGrad_i, srcGrad_j) > magnitude(dstGrad_i, dstGrad_j);
    drivingGrad_i = where(cond, srcGrad_i, dstGrad_i);
    drivingGrad_j = where(cond, srcGrad_j, dstGrad_j);
end

function [result] = magnitude(gx, gy)
    result = sqrt(gx.^2 + gy.^2);
end

function [result] = where(cond, x, y)
    result = zeros(size(cond));
    result(cond) = x(cond);
    result(~cond) = y(~cond);
end

function [translated] = get_translated_src(src, mask_src, mask_dst)
    [src_row,src_col,~] = find(mask_src == 1);
    [dst_row,dst_col,~] = find(mask_dst == 1);

    i_diff = dst_row(1,1)-src_row(1,1);
    j_diff = dst_col(1,1)-src_col(1,1);

    translated = imtranslate(src,[j_diff, i_diff]);
end