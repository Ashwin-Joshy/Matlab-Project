[im, map] = imread('D:\project\MAIN PROJECT\Matlab-Project\TestImages\Lena.tiff');
[j,k]=rgb2gray(I)
imshow(gray2rgb(j,k),[])
function [grayImage,colormap] = rgb2gray1(rgbImage)
[grayImage,colormap] = in2gray(rgbImage);
end
function rgbImage = gray2rgb(grayImage,colormap)
rgbImage = ind2rgb(grayImage, colormap);
end