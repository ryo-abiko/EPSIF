%{
Please cite the following paper if you use this code:

Abiko Ryo, and Masaaki Ikehara. 
"Fast Edge Preserving 2D Smoothing Filter Using Indicator Function." 
ICASSP 2019-2019 IEEE International Conference on Acoustics, 
Speech and Signal Processing (ICASSP). IEEE, 2019.

output = EPSIF(input,sigma_i,w_size,iter_max,mode,sigma_s,sigma_r,guide)

sigma_i  : parameter of indicator function (default:0.45)
w_size   : filter size (odd) (dafault: 9)
iter_max : number of iterations (default: 3)
mode     : filtering mode --> 'simple' Simple moving average (default)
                              'full'   bilateral filter
                              'guided' guided filtering
sigma_s  : spatial smoothing parameter (only for 'full')
sigma_r  : color range smoothing parameter (only for 'full')
guide    : guide image (only for 'guide')
%--------------------------------------------------------------------------
%}

clear all; close all;
addpath utilities/

image = im2double(imread('bird.png'));

smoothed_image = EPSIF(image);

figure(1)
imshow(smoothed_image)
imwrite(smoothed_image,'smoothed_image.png')

