function filter = filter_logGabor_derivative(imSize,degSize,freq,orientation,bw,f,orient0)
%function filter = filter_logGabor(imSize,degSize,freq,orientation,bw)
% This function creates a log-Gabor-derivative filter in frequency space. It takes an
% image size in px and degrees, a frequency in cyc/deg an
% orientation(0-2pi) and a bandwidth [freq,orientation] as input.
% It then computes the derivative with respect to the bandwidths bw(1:2)


global useGPU
if isempty(useGPU)
    useGPU = false;
end

if isempty(f) || isempty(orient0)
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
    filter = zeros([imSize,length(orientation),1,2],'gpuArray');
    pie = gpuArray(pi);
else
    filter = zeros([imSize,length(orientation),1,2]);
    pie = pi;
end

for iOrient = 1:length(orientation)
    orient = orient0-orientation(iOrient)+pie;
    orient = mod(orient,2*pie)-pie;
    filter(:,:,iOrient,1,1) = 2* (log2(f)-log2(freq)).^2./bw(1).^3.*exp(-1/2*((log2(f)-log2(freq)).^2./bw(1).^2+(orient.^2./bw(2)^2)));
    filter(ceil((size(f,1)+1)/2),ceil((size(f,2)+1)/2),iOrient,1,1) = 0; % no change at 0 frequency 
    filter(:,:,iOrient,1,2) = 2* (orient.^2./bw(2)^3) .* exp(-1/2*((log2(f)-log2(freq)).^2./bw(1).^2+(orient.^2./bw(2)^2)));
end