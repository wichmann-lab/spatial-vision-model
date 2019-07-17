function imageBitmap = image_disk(imSize,degSize,bgLum,lum,center,size)
%function image = image_disk(imSize,degSize,bgLum,lum,center,size)
% creates an image of a light disk lum brighter than the background at
% bgLum with center center and radius size in degrees.
% As usual it requires the size in pixels and in degrees.



x = linspace(0, degSize(2),imSize(2))-center(2);
y = linspace(0, degSize(1),imSize(1))-center(1);

r = bsxfun(@plus,x.^2,y'.^2);

imageBitmap = bgLum.*ones(imSize)+lum.* (r<size.^2);

