function [diff, noise, p,pGrad,noiseGrad,diffGrad] = compare_images(image1,image1Noise,image2,image2Noise,image1Grad,image1NoiseGrad,image2Grad,image2NoiseGrad)
% function [diff, noise, p] = compare_images(image1,image1_noise,image2,image2_noise,type)
% this function produces a comparison of two images.
% it returns the mean difference of pixels and its noise plus a probability
% of detection assuming a single draw with a fixed optimally placed
% criterion, i.e. one 2AFC trial.
%
% This can also calculate the gradient if the gradients of the images and
% their noise are also given.
% However this works only for the liklihood at the moment! This will
% produce an error if the gradients are not provided, or the wrong type is
% chosen.
% The gradient calculated here is actually the gradient of the overall
% signal to noise ratio!

global useGPU

if isempty(useGPU)
    useGPU = 0;
end

if ~exist('type','var') || isempty(type)
    type = 1;
end

if nargout >3
    assert(type ==1 && exist('image1Grad','var') && exist('image1NoiseGrad','var') && exist('image2Grad','var') && exist('image2NoiseGrad','var'),'To calculate a gradient you need to provide gradients from earlier processing!');
end

imDiff = image1-image2;
SNR  = imDiff(:)./sqrt(image1Noise(:)+image2Noise(:));
diff = sum(SNR.*SNR); % weigthed by signal to noise ratio, assuming independet pixels
noise = sum(SNR.*SNR);
if useGPU
    diff = gather(diff);
    noise = gather(noise);
end
p = normcdf(diff,0,sqrt(noise));

if nargout >3
    SNR = reshape(SNR,size(imDiff));
    %diffGrad = sum(sum(sum(sum(bsxfun(@times,bsxfun(@times,2.*imDiff./sqrt(image1Noise+image2Noise),(image1Grad-image2Grad))-bsxfun(@times,imDiff.^2./sqrt(image1Noise+image2Noise).^3,image1NoiseGrad+image2NoiseGrad),sign(imDiff)),1),2),3),4);
    diffGrad = bsxfun(@times,SNR,(image1Grad-image2Grad)-bsxfun(@times,imDiff./4./(image1Noise+image2Noise),(image1NoiseGrad+image2NoiseGrad)));
    diffGrad = 2*sum(sum(sum(sum(diffGrad,1),2),3),4);
    noiseGrad = sum(sum(sum(sum(bsxfun(@times,2.*imDiff,(image1Grad-image2Grad)),1),2),3),4);
    pGrad = 1./sqrt(noise).*diffGrad(:) - 0.5 * diff.*noise.^(-3/2).*noiseGrad(:);
    pGrad = normpdf(diff./sqrt(noise)).*pGrad;
end

if diff == 0
    p = .5;
end