% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

clc
clear all
close all

I=rgb2gray(importdata('A380small.jpg'));
imshow(I)

I2=imcrop(I,[25 100 50 50]);

imshow(I2)

y = xcorr2(double(I2),double(I));
figure; imagesc(y)