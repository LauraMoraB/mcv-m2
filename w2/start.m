%% Generate masks
clearvars;
dst = imread('images/Atardecer.jpg');
src = imread('images/bird_800x450.jpg');

[mask_src] = getMask(src);
[~, mask_dst] = alignSource(src, mask_src, dst);

mask_src=logical(mask_src);
mask_dst=logical(mask_dst);
src=double(src);
dst=double(dst);

[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    [drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    %[drivingGrad_i, drivingGrad_j] = mixGradients(src(:,:,nC), dst(:,:,nC), mask_src, mask_dst);
    
    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);
     
    driving_on_dst = zeros(size(dst(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;
    [dst1(:,:,nC),t] = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst, param);
    
end
imshowpair(dst/256, dst1/256, 'montage')
%imwrite(dst1, 'results/dst_src_name.png')

%% LENA and GIRL
dst = double(imread('images/lena.png'));
src = double(imread('images/girl.png')); % flipped girl, because of the eyes
[ni,nj, nChannels]=size(dst);

param.hi=1;
param.hj=1;


%Eyes
mask_src=logical(imread('images/mask_src_eyes.png'));
mask_dst=logical(imread('images/mask_dst_eyes.png'));

for nC = 1: nChannels
    %TO DO: Compute Gradient to solve the Poisson Equation
    
    %[drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    [drivingGrad_i, drivingGrad_j] = mixGradients(src(:,:,nC), ...
                                                dst(:,:,nC), ...
                                                mask_src, mask_dst);

    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);

    driving_on_dst = zeros(size(src(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));

    param.driving = driving_on_dst;

    [dst1(:,:,nC), t] = sol_Poisson_Equation_Axb(dst(:,:,nC), mask_dst, param);
end

%Mouth
%masks to exchange: Mouth
mask_src=logical(imread('images/mask_src_mouth.png'));
mask_dst=logical(imread('images/mask_dst_mouth.png'));

for nC = 1: nChannels
    
    %TO DO: COMPLETE the ??
    %[drivingGrad_i, drivingGrad_j] = importGradients(src(:,:,nC));
    [drivingGrad_i, drivingGrad_j] = mixGradients(src(:,:,nC), dst(:,:,nC), mask_src, mask_dst);

    driving_on_src = divergence(drivingGrad_i, drivingGrad_j);
    
    driving_on_dst = zeros(size(src(:,:,1)));
    driving_on_dst(mask_dst(:)) = driving_on_src(mask_src(:));
    
    param.driving = driving_on_dst;
    [dst1(:,:,nC), t] = sol_Poisson_Equation_Axb(dst1(:,:,nC), mask_dst, param);
end

figure(3);
imshowpair(dst/256, dst1/256, 'montage')
%imwrite(dst1, 'results/lena_girl.png')

%% Auxiliary functions
function [drivingGrad_i, drivingGrad_j] = importGradients(src)
    [srcGrad_i, srcGrad_j] = gradient(src);
    drivingGrad_i = srcGrad_i;
    drivingGrad_j = srcGrad_j;
end

function [drivingGrad_i, drivingGrad_j] = mixGradients(src, dst, ...
                                                        mask_src, mask_dst)
    % Compute gradient by applying mixing gradient technique
    [srcGrad_i, srcGrad_j] = gradient(src);
    [dstGrad_i, dstGrad_j] = gradient(dst);
    srcGrad_i = srcGrad_i(mask_src(:));
    srcGrad_j = srcGrad_j(mask_src(:));
    dstGrad_i = dstGrad_i(mask_dst(:));
    dstGrad_j = dstGrad_j(mask_dst(:));
    cond = magnitude(srcGrad_i, srcGrad_j) > magnitude(dstGrad_i, dstGrad_j);
    drivingGrad_i = zeros(size(src));
    drivingGrad_j = zeros(size(src));
    drivingGrad_i(mask_src(:)) = where(cond, srcGrad_i, dstGrad_i);
    drivingGrad_j(mask_src(:)) = where(cond, srcGrad_j, dstGrad_j);
end

function [result] = magnitude(gx, gy)
    result = sqrt(gx.^2 + gy.^2);
end

function [result] = where(cond, x, y)
    % apply condition. Keep gradient with higher magnitude
    result = zeros(size(cond));
    result(cond) = x(cond);
    result(~cond) = y(~cond);
end