function imageBitmapPad = pad_image(imageBitmap, type)
% function pad_image(image, type)
% this function pads an image to the next power of two. There are different
% types of doing so. 
% 0 = pad with mean value;
% 1 = padding with border values
% 2 = pad with zeros

global useGPU
if isempty(useGPU)
    useGPU = false;
end

imSize = 2.^(ceil(log2(max(size(imageBitmap)))));
if useGPU
    imageBitmapPad = mean(imageBitmap(:)).*gpuArray.ones(imSize); 
else
    imageBitmapPad = mean(imageBitmap(:)).*ones(imSize); 
end
xID = floor((imSize-size(imageBitmap,2))/2);
xID = [xID+1,xID+size(imageBitmap,2)];
yID = floor((imSize-size(imageBitmap,1))/2);
yID = [yID+1,yID+size(imageBitmap,1)];
switch type
    case 0
        % pad with mean value;
        imageBitmapPad(yID(1):yID(2),xID(1):xID(2)) = imageBitmap;
    case 1
        imageBitmapPad(yID(1):yID(2),xID(1):xID(2)) = imageBitmap;
        % padding with border values
        imageBitmapPad(1:yID(1),:)  = repmat(imageBitmapPad(yID(1),:),[yID(1),1]);
        imageBitmapPad(:,1:xID(1))  = repmat(imageBitmapPad(:,xID(1)),[1,xID(1)]);
        imageBitmapPad(yID(2):end,:)= repmat(imageBitmapPad(yID(2),:),[imSize-yID(2)+1,1]);
        imageBitmapPad(:,xID(2):end)= repmat(imageBitmapPad(:,xID(2)),[1,imSize-xID(2)+1]);
    case 2
        % pad with zeros
        imageBitmapPad  = zeros(imSize);
        imageBitmapPad(yID(1):yID(2),xID(1):xID(2)) = imageBitmap;
    otherwise
        error('given padding not recognized');
end