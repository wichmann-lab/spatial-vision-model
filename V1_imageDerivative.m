function imgD = V1_imageDerivative(outD,imageBitmap,degSize,pars)
%function imgD = V1_imageDerivative(outD,imageBitmap,degSize,pars)
%
% This is for calculating the derivative of the detectablitly in image
% space.
% It takes outD as input which should be the derivative d(d')/d(dout_i)


global useGPU

if isempty(useGPU)
    useGPU = false;
end

assert(length(pars)>=12,'You need to provide a std for orientation normalization pool')

% pars are fixed in the model files and should not default here
CNaka = pars(1);    % Naka rushton constant
ExNaka = pars(2)+pars(3);   % Naka rushton Exponent top
ExNakaNorm = pars(2); % Naka rushton Exponent below
bw(1) = pars(4);    % bandwidth in frequency (std of log-Gabor in octaves)
bw(2) = pars(5);    % bandwidth in orientation (std log-Gabor in radiants)
nFreq = pars(6);    % number of frequency bands
nOrient = pars(7);  % number of orientations
poolSize = pars(8); % meaning depends on type, spatial size of normalization pool
poolbw = pars(9);   % bandwidth of the normalization pool (Octaves in frequency)
minF   = pars(10);  % lowest frequency band
maxF   = pars(11);  % highest frequency band
poolO = pars(12);

if useGPU
    imageBitmap = gpuArray(imageBitmap);
end

[out,frequencies,~,~,filters] = decomp_Gabor(imageBitmap,degSize,[minF,maxF],nFreq,nOrient,bw);
ao = abs(out);

imSize = size(imageBitmap);

%% MEAN NORMALIZATION -> type 6 
lao = log(ao);
%lao0 = lao;
%lao0(ao==0) = 0;
%normalizer1 = abs(out).^ExNakaNorm;
normalizer1 = exp(lao.*ExNakaNorm);
normalizer = mean(mean(normalizer1,1),2);
fdiff = linspace(log2(minF)-log2(maxF),log2(maxF)-log2(minF),2*length(frequencies)-1);
if useGPU
    fdiff = gpuArray(fdiff);
end
gaussF = exp(-fdiff.^2./poolbw.^2./2);
gaussFN = gaussF./sum(gaussF(:));
gaussFN = reshape(gaussFN,1,1,1,[]);
o = ((0:(size(out,3)-1))-floor(size(out,3)/2))*pi/size(out,3);
if useGPU
    o = gpuArray(o);
end
if ~isfinite(poolO)
    gaussO = ones(size(o));
elseif poolO == 0
    gaussO = (o == 0);
else
    gaussO  = exp(-o.^2./poolO.^2./2);
end
gaussON = gaussO./sum(gaussO);
gaussONF = fft(ifftshift(gaussON));
%gaussO = reshape(gaussO,1,1,[]);
gaussONF = gaussONF(:);
% fourier space filtering to enable "wrap around"
normalizerF= fft(squeeze(normalizer));
normalizerF = bsxfun(@times,normalizerF,gaussONF);
normalizer = ifft(normalizerF,'symmetric');
normalizer = reshape(normalizer,1,1,nOrient,nFreq);
normalizerPad = padarray(normalizer,[0,0,0,size(normalizer,4)-1]);
normalizer = convn(normalizerPad,gaussFN,'valid');
normalizerFin = (normalizer+CNaka.^ExNakaNorm);
%outNormalized = bsxfun(@rdivide,abs(out).^ExNaka,normalizerFin);
% only necessary for testing
%outNormalized = bsxfun(@rdivide,exp(lao.*ExNaka),normalizerFin);

%% laoD
laoD = bsxfun(@rdivide,ExNaka.*exp(lao.*(ExNaka)),normalizerFin);

normalizerFinD = -sum(sum(bsxfun(@rdivide,outD.*exp(lao.*ExNaka),normalizerFin.^2),1),2);

normalizerFinDPad = padarray(normalizerFinD,[0,0,0,size(normalizerFinD,4)-1]);
dnormalizer = convn(normalizerFinDPad,gaussFN,'valid'); % is symmetric -> no adjustment needed

dnormalizerF = fft(squeeze(dnormalizer));
dnormalizerF = bsxfun(@times,dnormalizerF,fft(ifftshift(flip(gaussON(:)))));
dnormalizer  = ifft(dnormalizerF,'symmetric');
dnormalizer = repmat(reshape(dnormalizer./prod(imSize),1,1,nOrient,nFreq),imSize);
dnormalizer1 = dnormalizer.*ExNakaNorm.*exp(lao.*(ExNakaNorm));

% convert Filters
if useGPU
    pyrR = gpuArray.nan(size(outD));
    pyrI = gpuArray.nan(size(outD));
else
    pyrR = nan(size(outD));
    pyrI = nan(size(outD));
end
for iFilter = 1:nFreq
    filtInvR = conj(fft2(real(ifft2(filters{iFilter}))));
    filtInvI = conj(fft2(imag(ifft2(filters{iFilter}))));
    pyrR(:,:,:,iFilter) = ifft2(filtInvR.*fft2(real(dout(:,:,:,iFilter))),'symmetric');
    pyrI(:,:,:,iFilter) = ifft2(filtInvI.*fft2(imag(dout(:,:,:,iFilter))),'symmetric');
end
 
imgD = sum(sum(pyrR+pyrI,4),3);

