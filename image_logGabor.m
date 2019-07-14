function image = image_logGabor(imSize,degSize,freq,contrast,center,orientation,L,bw,phase,antiAliasing)
% function image = image_logGabor(imSize,degSize,freq,contrast,center,orientation,L,bw,phase)
% this function creates a Gabor patch with given orientation, frequency and
% size (rad and deg of visual angle) in the center of the image
%
% 0 phase is the odd positive phase
%
% numerical normalization to 1 at peak of the Gabor


global useGPU
if isempty(useGPU)
    useGPU = false;
end

if ~exist('antiAliasing','var') || isempty(antiAliasing)
    antiAliasing = 4;
end
imSize = antiAliasing*imSize;


if ~exist('f','var') || isempty(f) || ~exist('orient0','var') || isempty(orient0)
    if useGPU
        x = gpuArray.linspace(-imSize(2)/2/degSize(2),(imSize(2)-2)/2/degSize(2),imSize(2));
        y = gpuArray.linspace(-imSize(1)/2/degSize(1),(imSize(1)-2)/2/degSize(1),imSize(1));
    else
        x = 1:imSize(2);
        y = 1:imSize(1);
        x = (x-ceil(mean(x)))./degSize(2);
        y = (y-ceil(mean(y)))./degSize(1);
    end
    f = sqrt(bsxfun(@plus,x.^2,y'.^2));
    orient0 = bsxfun(@(x,y) atan2(y,x),x,y');
end

if useGPU
    pie = gpuArray(pi);
    orient = orient0-orientation+pie;
    %orient(orient>pie)  = orient(orient>pie)-2*pie;
    %orient(orient<-pie) = orient(orient<-pie)+2*pie;
    orient = mod(orient,2*pie)-pie;
    
    grating = 2* exp(-(log2(f)-log2(freq)).^2./bw(1).^2-(orient.^2./bw(2)^2)); % bw is two standard deviations
else
    orient = orient0-orientation;
    orient(orient>pi) = orient(orient>pi)-2*pi;
    orient(orient<-pi) = orient(orient<-pi)+2*pi;
    
    grating = 2* exp(-(log2(f)-log2(freq)).^2./bw(1).^2-(orient.^2./bw(2)^2)); % bw is two standard deviations
end

complexGabor = fftshift(ifft2(ifftshift(grating)));
Gabor = 1./max(abs(complexGabor(:))).*(cos(phase).*real(complexGabor)+sin(phase).*imag(complexGabor));

center0 = degSize/2;
shiftDeg = center-center0;
shiftPix = imSize./degSize.*shiftDeg;
if shiftPix(1) > 0
    Gabor((round(shiftPix(1))+1):end,:) = Gabor(1:(end-(round(shiftPix(1)))),:);
elseif shiftPix(1)<0
    Gabor(1:(end+(round(shiftPix(1)))),:) = Gabor((round(-shiftPix(1))+1):end,:);
end
if shiftPix(2) > 0
    Gabor(:,(round(shiftPix(2))+1):end) = Gabor(:,1:(end-(round(shiftPix(2)))));
elseif shiftPix(2)<0
    Gabor(:,1:(end+(round(shiftPix(2))))) = Gabor(:,(round(-shiftPix(2))+1):end);
end

image = L + contrast.*Gabor.*L;

image = imresize(image,1/antiAliasing);