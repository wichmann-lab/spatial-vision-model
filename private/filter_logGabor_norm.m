function [filter,NormFilterReal,NormFilterImag] = filter_logGabor_norm(imSize,degSize,freq,orientation,bw,f,orient0,NormSize,NormBw)

% Be aware that this is with wrap around... can go wrong for large
% Normalizer sizes!
% THIS DOES NOT WORK IMEADIATELY! FILTERS CANCEL EACH OTHER!


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
    filter = zeros([imSize,length(orientation),length(freq)],'gpuArray');
    orient = zeros([imSize],'gpuArray');
else
    filter = zeros([imSize,length(orientation),length(freq)]);
end

if useGPU
    pie = gpuArray(pi);
    for iFreq = 1:length(freq)
        for iOrient = 1:length(orientation)
            orient = orient0-orientation(iOrient)+pie;
            %orient(orient>pie)  = orient(orient>pie)-2*pie;
            %orient(orient<-pie) = orient(orient<-pie)+2*pie;
            orient = mod(orient,2*pie)-pie;
            
            filter(:,:,iOrient,iFreq) = 2* exp(-0.5*((log2(f)-log2(freq(iFreq))).^2./bw(1).^2+(orient.^2./bw(2)^2)));
        end
    end
else
    for iFreq = 1:length(freq)
        for iOrient = 1:length(orientation)
            orient = orient0-orientation(iOrient)+pi;
            %         orient(orient>pi) = orient(orient>pi)-2*pi;
            %         orient(orient<-pi) = orient(orient<-pi)+2*pi;
            orient = mod(orient,2*pi)-pi;
            
            filter(:,:,iOrient,iFreq) = 2* exp(-0.5*((log2(f)-log2(freq(iFreq))).^2./bw(1).^2+(orient.^2./bw(2)^2)));
        end
    end
end

%% set leftover high frequency to 0 for even filter size
if ~mod(size(filter,1),2)
    filter(1,:,:,:) = 0;
end
if ~mod(size(filter,2),2)
    filter(:,1,:,:) = 0;
end


%% calculate filters for Normalizer

if nargout >= 2
    % calculate real and imaginary filter
    if all(~mod(size(filter(:,:,1)),2)) % all even size
        realFilter = filter;
        realFilter(2:end,2:end,:) = 0.5*filter(2:end,2:end,:)+0.5*flip(flip(filter(2:end,2:end,:),1),2);
        imagFilter = filter;
        imagFilter(2:end,2:end,:) = 1i.*(-0.5*filter(2:end,2:end,:)+0.5*flip(flip(filter(2:end,2:end,:),1),2));
    end
    % calculate spatial Pooling
    if NormSize(1)==NormSize(2) % i.e. We can reuse f
        SpacePool = exp(-(f.*NormSize(1)).^2.*pi.^2.*2);
    else
        f = sqrt(bsxfun(@plus,(NormSize(2).*x).^2,(NormSize(1).*y)'.^2));
        SpacePool = exp(-f.^2.*pi.^2.*2);
    end
    realFilter = bsxfun(@times,realFilter,SpacePool);
    imagFilter = bsxfun(@times,imagFilter,SpacePool);
    % calculate orientation pooling
    for iOrient = 1:length(orientation)
        weights = mod((orientation-orientation(iOrient)+pi/2),pi)-pi/2;
        weights = exp(-weights.^2./2./NormBw(2).^2);
        weights = reshape(weights,1,1,[]);
        realFilter(:,:,iOrient,:)= sum(bsxfun(@times,weights,realFilter),3);
        imagFilter(:,:,iOrient,:)= sum(bsxfun(@times,weights,imagFilter),3);
    end
    % calculate frequency pooling
    for iFreq = 1:length(freq)
        weights = log2(freq)-log2(freq(iFreq));
        weights = exp(-weights.^2./2./NormBw(1).^2);
        weights = reshape(weights,1,1,1,[]);
        realFilter(:,:,:,iFreq)= sum(bsxfun(@times,weights,realFilter),4);
        imagFilter(:,:,:,iFreq)= sum(bsxfun(@times,weights,imagFilter),4);
    end
end