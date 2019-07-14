function filter = filter_logGabor(imSize,degSize,freq,orientation,bw,f,orient0)
%function filter = filter_logGabor(imSize,degSize,freq,orientation,bw)
% This function creates a log-Gabor filter in frequency space. It takes an
% image size in px and degrees, a frequency in cyc/deg an
% orientation(0-2pi) and a bandwidth [freq,orientation] as input.


global useGPU
if isempty(useGPU)
    useGPU = false;
end

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
    filter = zeros([imSize,length(orientation)],'gpuArray');
    %orient = zeros([imSize],'gpuArray');
else
    filter = zeros([imSize,length(orientation)]);
end

if useGPU
    pie = gpuArray(pi);
    for iOrient = 1:length(orientation)
        orient = orient0-orientation(iOrient)+pie;
        %orient(orient>pie)  = orient(orient>pie)-2*pie;
        %orient(orient<-pie) = orient(orient<-pie)+2*pie;
        orient = mod(orient,2*pie)-pie;
        
        filter(:,:,iOrient) = 2* exp(-1/2*((log2(f)-log2(freq)).^2./bw(1).^2+(orient.^2./bw(2)^2))); 
    end
else
    for iOrient = 1:length(orientation)
        orient = orient0-orientation(iOrient)+pi;
        %orient(orient>pi) = orient(orient>pi)-2*pi;
        %orient(orient<-pi) = orient(orient<-pi)+2*pi;
        orient = mod(orient,2*pi)-pi;
        
        filter(:,:,iOrient) = 2* exp(-1/2*((log2(f)-log2(freq)).^2./bw(1).^2+(orient.^2./bw(2)^2))); 
    end
end
%filter = fftshift(filter,1);
%filter = fftshift(filter,2);