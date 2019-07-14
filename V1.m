function [outNormalized,outGrad,out] = V1(imageBitmap,degSize,V1Mode,pars,Gradidx)
%function [out,outGrad] = V1(imageBitmap,degSize,type,pars,Gradidx)
%
% this runs a log-Gabor pyramid decompostion and computes a normalization
% of the different channels.
%
% Gradient is only implemented for number 4+6, i.e. the constant size
% normalization pool variants
%
% There is a severe bug in the current gpuArray.ifft implementation for
% symmetric. I cured this for number 4 and 6 but all other variants of the model
% remain unchecked: Be freakishly careful with this....
%
% Derivatives can be restricted to only some parameters using the Gradidx
% argument. It should then contain the parameters you want the derivatives
% to.


global useGPU

if isempty(useGPU)
    useGPU = false;
end

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

if length(pars)>=12
    poolO = pars(12);
end


if nargout >= 2
    if ~exist('Gradidx','var') 
        Gradidx = [1,2,3,4,5,9,12];
    end
    GradOrder = nan(size(pars));
    for iPar = 1:length(pars)
        if ismember(iPar,Gradidx)
            GradOrder(iPar) = find(Gradidx==iPar);
        end
    end
    
    if useGPU
        outGrad = gpuArray.zeros([size(imageBitmap),nOrient,nFreq,length(Gradidx)]);
    else
        outGrad = zeros([size(imageBitmap),nOrient,nFreq,length(Gradidx)]);
    end
end


%% first independent of the type we calculate the frequency decomposition
% into log-Gabor frequency and orientation bands
if useGPU
    imageBitmap = gpuArray(imageBitmap);
end

if nargout >=2 && any(~isnan(GradOrder(4:5)))
    [out,frequencies,pyrGradBW] = decomp_Gabor(imageBitmap,degSize,[minF,maxF],nFreq,nOrient,bw);
    ao = abs(out);
    if all(~isnan(GradOrder(4:5)))
        outGrad(:,:,:,:,GradOrder(4:5)) = bsxfun(@rdivide,bsxfun(@times,real(out),real(pyrGradBW))+bsxfun(@times,imag(out),imag(pyrGradBW)),max(ao,eps));
    elseif ~isnan(GradOrder(4))
        outGrad(:,:,:,:,GradOrder(4)) = bsxfun(@rdivide,bsxfun(@times,real(out),real(pyrGradBW(:,:,:,:,1)))+bsxfun(@times,imag(out),imag(pyrGradBW(:,:,:,:,1))),max(ao,eps));
    elseif ~isnan(GradOrder(5))
        outGrad(:,:,:,:,GradOrder(5)) = bsxfun(@rdivide,bsxfun(@times,real(out),real(pyrGradBW(:,:,:,:,2)))+bsxfun(@times,imag(out),imag(pyrGradBW(:,:,:,:,2))),max(ao,eps));
    end
else
    [out,frequencies] = decomp_Gabor(imageBitmap,degSize,[minF,maxF],nFreq,nOrient,bw);
    ao = abs(out);
end

imSize = size(imageBitmap);
x = 1:imSize(2);
y = 1:imSize(1);
x = (x-ceil(mean(x))).*degSize(2)./imSize(2);
y = (y-ceil(mean(y))).*degSize(1)./imSize(1);

if nargout>= 2
    r = bsxfun(@plus,x.^2,y'.^2);
end


switch V1Mode
    case 1 % spatially pool whole image (Probably only sensible for foveal model)
        % For obvious reasons this ignores the size of the pool in space
        assert(length(pars)>=12,'You need to provide a std for orientation normalization pool')
        lao = log(ao);
        lao0 = lao;
        lao0(ao==0) = 0;
        %normalizer1 = abs(out).^ExNakaNorm;
        normalizer1 = exp(lao.*ExNakaNorm);
        normalizer = mean(mean(normalizer1,1),2);
        fdiff = linspace(log2(minF)-log2(maxF),log2(maxF)-log2(minF),2*length(frequencies)-1);
        gaussF = exp(-fdiff.^2./poolbw.^2./2);
        gaussFN = gaussF./sum(gaussF(:));
        gaussFN = reshape(gaussFN,1,1,1,[]);
        o = ((0:(size(out,3)-1))-floor(size(out,3)/2))*pi/size(out,3);
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
        if nargout>=2
            if ~isnan(GradOrder(9))
                gaussFD = exp(-fdiff.^2./poolbw.^2./2).*fdiff.^2./poolbw.^3;
                gaussFD = (gaussFD(:)- gaussFN(:).*sum(gaussFD(:)))./sum(gaussF(:));
                gaussFD = reshape(gaussFD,1,1,1,[]);
            end
            if ~isnan(GradOrder(12))
                gaussOD = exp(-o.^2./poolO.^2./2).*o.^2./poolO.^3;
                gaussOD = (gaussOD(:)- gaussON(:).*sum(gaussOD(:)))./sum(gaussO(:));
                gaussODF = fft(ifftshift(gaussOD));
                gaussODF = gaussODF(:);
                BdO     = bsxfun(@times,normalizerF,gaussODF);
                BdO     = ifft(BdO,'symmetric');
            end
        end
        normalizerF = bsxfun(@times,normalizerF,gaussONF);
        normalizer = ifft(normalizerF,'symmetric');
        normalizer = reshape(normalizer,1,1,nOrient,nFreq);
        normalizerPad = padarray(normalizer,[0,0,0,size(normalizer,4)-1]);
        if nargout>=2
            if ~isnan(GradOrder(12))
                BdO = reshape(BdO,1,1,nOrient,nFreq);
                BdOPad = padarray(BdO,[0,0,0,size(BdO,4)-1]);
                BdO = convn(BdOPad,gaussFN,'valid');
            end
            if ~isnan(GradOrder(9))
                BdF     = convn(normalizerPad,gaussFD,'valid');
            end
        end
        normalizer = convn(normalizerPad,gaussFN,'valid');
        normalizerFin = (normalizer+CNaka.^ExNakaNorm);
        %outNormalized = bsxfun(@rdivide,abs(out).^ExNaka,normalizerFin);
        outNormalized = bsxfun(@rdivide,exp(lao.*ExNaka),normalizerFin);
        if nargout >= 2 % calculate gradient
            if ~isnan(GradOrder(2))
                logNormalizer = lao0.*normalizer1;
                logNormalizer = mean(mean(logNormalizer,1),2);
                logNormalizerF= fft(squeeze(logNormalizer));
                logNormalizerF = bsxfun(@times,logNormalizerF,gaussONF);
                logNormalizer = ifft(logNormalizerF,'symmetric');
                logNormalizer = reshape(logNormalizer,1,1,nOrient,nFreq);
                logNormalizerPad = padarray(logNormalizer,[0,0,0,size(logNormalizer,4)-1]);
                logNormalizer = convn(logNormalizerPad,gaussFN,'valid');
            end
            
            if ~isnan(GradOrder(4))
                Bdf = exp(lao.*(ExNakaNorm-1)).*outGrad(:,:,:,:,GradOrder(4));
                Bdf = ExNakaNorm*mean(mean(Bdf,1),2);
                BdfF= fft(squeeze(Bdf));
                BdfF = bsxfun(@times,BdfF,gaussONF);
                Bdf = ifft(BdfF,'symmetric');
                Bdf = reshape(Bdf,1,1,nOrient,nFreq);
                BdAPad = padarray(Bdf,[0,0,0,size(Bdf,4)-1]);
                Bdf = convn(BdAPad,gaussFN,'valid');
            end
            
            if ~isnan(GradOrder(5))
                Bdphi = exp(lao.*(ExNakaNorm-1)).*outGrad(:,:,:,:,GradOrder(5));
                Bdphi = ExNakaNorm*mean(mean(Bdphi,1),2);
                BdphiF= fft(squeeze(Bdphi));
                BdphiF = bsxfun(@times,BdphiF,gaussONF);
                Bdphi = ifft(BdphiF,'symmetric');
                Bdphi = reshape(Bdphi,1,1,nOrient,nFreq);
                BdAPad = padarray(Bdphi,[0,0,0,size(Bdphi,4)-1]);
                Bdphi = convn(BdAPad,gaussFN,'valid');
            end
            
            if any(~isnan(GradOrder([4,5,9,12])))
                RdB = bsxfun(@rdivide,-outNormalized,normalizerFin);
            end
            if any(~isnan(GradOrder([4,5])))
                RdA = ExNaka.* bsxfun(@rdivide,exp(lao.*(ExNaka-1)),normalizerFin);
            end
            
            if ~isnan(GradOrder(1))
                outGrad(:,:,:,:,GradOrder(1)) = bsxfun(@rdivide,-outNormalized.*ExNakaNorm.*(CNaka.^(ExNakaNorm-1)),normalizerFin);
            end
            if ~isnan(GradOrder(2))
                outGrad(:,:,:,:,GradOrder(2)) = outNormalized.*(bsxfun(@minus,lao0, bsxfun(@rdivide,logNormalizer+log(CNaka).*CNaka.^ExNakaNorm,normalizerFin)));
            end
            if ~isnan(GradOrder(3))
                outGrad(:,:,:,:,GradOrder(3)) = lao0.*outNormalized;
            end
            if ~isnan(GradOrder(4))
                outGrad(:,:,:,:,GradOrder(4)) = RdA.*outGrad(:,:,:,:,GradOrder(4))+bsxfun(@times,RdB,Bdf);
            end
            if ~isnan(GradOrder(5))
                outGrad(:,:,:,:,GradOrder(5)) = RdA.*outGrad(:,:,:,:,GradOrder(5))+bsxfun(@times,RdB,Bdphi);
            end
            if ~isnan(GradOrder(9))
                outGrad(:,:,:,:,GradOrder(9))  = bsxfun(@times,RdB,BdF);
            end
            if ~isnan(GradOrder(12))
                outGrad(:,:,:,:,GradOrder(12)) = bsxfun(@times,RdB,BdO);
            end
        end
    case 2 % local normalization -> only activations at this pixel count
        assert(length(pars)>=12,'You need to provide a std for orientation normalization pool')
        
        lao = log(ao);
        lao0 = lao;
        lao0(ao==0) = 0;
        normalizer1 = exp(lao.*ExNakaNorm);
        
        normalizer = normalizer1;
                
        fdiff   = linspace(log2(minF)-log2(maxF),log2(maxF)-log2(minF),2*length(frequencies)-1);
        gaussF  = exp(-fdiff.^2./poolbw.^2./2);
        gaussFN = gaussF./sum(gaussF(:));
        gaussFN = reshape(gaussFN,1,1,1,[]);
        o = ((0:(size(out,3)-1))-floor(size(out,3)/2))*pi/size(out,3);
        if ~isfinite(poolO)
            gaussO = ones(size(o));
        elseif poolO == 0
            gaussO = (o == 0);
        else
            gaussO  = exp(-o.^2./poolO.^2./2);
        end
        gaussON = gaussO./sum(gaussO);
        gaussONF= fft(ifftshift(gaussON));
        gaussONF= reshape(gaussONF,1,1,[]);
        
        if useGPU
            gaussONF = gpuArray(gaussONF);
            gaussFN = gpuArray(gaussFN);
        end
        % fourier space filtering to enable "wrap around"
        normalizerF = fft(normalizer,[],3);
        if nargout >= 2
            if ~isnan(GradOrder(9))
                gaussFD = exp(-fdiff.^2./poolbw.^2./2).*fdiff.^2./poolbw.^3;
                gaussFD = (gaussFD(:)- gaussFN(:).*sum(gaussFD(:)))./sum(gaussF(:));
                gaussFD = reshape(gaussFD,1,1,1,[]);
                if useGPU
                    gaussFD = gpuArray(gaussFD);
                end
            end
            if ~isnan(GradOrder(12))
                gaussOD = exp(-o.^2./poolO.^2./2).*o.^2./poolO.^3;
                gaussOD = (gaussOD(:)- gaussON(:).*sum(gaussOD(:)))./sum(gaussO(:));
                gaussODF= fft(ifftshift(gaussOD));
                gaussODF= reshape(gaussODF,1,1,[]);
                BdO     = bsxfun(@times,normalizerF,gaussODF);
                BdO     = ifft(BdO,[],3,'symmetric');
            end
            if ~isnan(GradOrder(8))
                BdS     = bsxfun(@times,fft(BdS,[],3),gaussONF);
                BdS     = ifft(BdS,[],3,'symmetric');
            end
        end
        normalizerF = bsxfun(@times,normalizerF,gaussONF);
        normalizer  = ifft(normalizerF,[],3,'symmetric');
        clear normalizerF
        normalizerPad = padarray(normalizer,[0,0,0,size(normalizer,4)-1]);
        if nargout>=2
            if ~isnan(GradOrder(12))
                BdOPad = padarray(BdO,[0,0,0,size(BdO,4)-1]);
                BdO = convn(BdOPad,gaussFN,'valid');
                clear BdOPad
            end
            if ~isnan(GradOrder(9))
                BdF = convn(normalizerPad,gaussFD,'valid');
            end
            if ~isnan(GradOrder(8))
                BdSPad = padarray(BdS,[0,0,0,size(BdS,4)-1]);
                BdS = convn(BdSPad,gaussFN,'valid');
                clear BdSPad
            end
        end
        normalizer = convn(normalizerPad,gaussFN,'valid');
        clear normalizerPad
        normalizerFin = (normalizer+CNaka.^ExNakaNorm);
        outNormalized = bsxfun(@rdivide,exp(lao.*ExNaka),normalizerFin);
        if nargout >= 2 % calculate gradient
            if ~isnan(GradOrder(2))
                logNormalizer = lao0.*normalizer1;
                logNormalizerF= fft(logNormalizer,[],3);
                logNormalizerF = bsxfun(@times,logNormalizerF,gaussONF);
                logNormalizer = ifft(logNormalizerF,[],3,'symmetric');
                logNormalizerPad = padarray(logNormalizer,[0,0,0,size(logNormalizer,4)-1]);
                logNormalizer = convn(logNormalizerPad,gaussFN,'valid');
            end
            
            if ~isnan(GradOrder(4))
                Bdf   = ExNakaNorm.*exp(lao.*(ExNakaNorm-1)).*outGrad(:,:,:,:,GradOrder(4));
                BdfF= fft(Bdf,[],3);
                BdfF = bsxfun(@times,BdfF,gaussONF);
                Bdf = ifft(BdfF,[],3,'symmetric');
                BdAPad = padarray(Bdf,[0,0,0,size(Bdf,4)-1]);
                Bdf = convn(BdAPad,gaussFN,'valid');
            end
            
            if ~isnan(GradOrder(5))
                Bdphi = ExNakaNorm.*exp(lao.*(ExNakaNorm-1)).*outGrad(:,:,:,:,GradOrder(5));
                BdphiF= fft(Bdphi,[],3);
                BdphiF= bsxfun(@times,BdphiF,gaussONF);
                Bdphi = ifft(BdphiF,[],3,'symmetric');
                BdAPad= padarray(Bdphi,[0,0,0,size(Bdphi,4)-1]);
                Bdphi = convn(BdAPad,gaussFN,'valid');
            end
            
            if any(~isnan(GradOrder([4,5,8,9,12])))
                RdB = bsxfun(@rdivide,-outNormalized,normalizerFin);
            end
            if any(~isnan(GradOrder([4,5])))
                RdA = ExNaka.* bsxfun(@rdivide,exp(lao.*(ExNaka-1)),normalizerFin);
            end
            
            if ~isnan(GradOrder(1))
                outGrad(:,:,:,:,GradOrder(1)) = bsxfun(@rdivide,-outNormalized.*ExNakaNorm.*(CNaka.^(ExNakaNorm-1)),normalizerFin);
            end
            if ~isnan(GradOrder(2))
                outGrad(:,:,:,:,GradOrder(2)) = outNormalized.*(bsxfun(@minus,lao0, bsxfun(@rdivide,logNormalizer+log(CNaka).*CNaka.^ExNakaNorm,normalizerFin)));
            end
            if ~isnan(GradOrder(3))
                outGrad(:,:,:,:,GradOrder(3)) = lao0.*outNormalized;
            end
            if ~isnan(GradOrder(4))
                outGrad(:,:,:,:,GradOrder(4)) = RdA.*outGrad(:,:,:,:,GradOrder(4))+bsxfun(@times,RdB,Bdf);
            end
            if ~isnan(GradOrder(5))
                outGrad(:,:,:,:,GradOrder(5)) = RdA.*outGrad(:,:,:,:,GradOrder(5))+bsxfun(@times,RdB,Bdphi);
            end
            if ~isnan(GradOrder(8))
                outGrad(:,:,:,:,GradOrder(8)) = bsxfun(@times,RdB,BdS);
            end
            if ~isnan(GradOrder(9))
                outGrad(:,:,:,:,GradOrder(9)) = bsxfun(@times,RdB,BdF);
            end
            if ~isnan(GradOrder(12))
                outGrad(:,:,:,:,GradOrder(12))= bsxfun(@times,RdB,BdO);
            end
        end
end
