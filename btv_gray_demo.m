%% Gray image experiment
clc; clear;
close all;
seed = 2025;
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
%% Path Information
addpath(genpath('MyCode'));
addpath(genpath('quality asses'));
%% Load Data
x_ori = double(imread(strcat(pwd,'\test006.png')));
x_ori = mat2gray(x_ori);
imshow(x_ori,[])
[n1,n2] = size(x_ori); N = [n1,n2];
sr    = 0.1;
%% Sampling
Phi = PermuteWHT_partitioned(n1*n2,1,sr);
b   = Phi*x_ori(:);
%% BTV for image recovery
opts.N = [n1,n2];
opts.lambda      = 1e-3;
opts.block_size  = 4;
xrec_btv  = btv_admm_gray(b,Phi,opts);
[psnr, ssim] = msqia(x_ori, xrec_btv)

% The result is [0.0714, 22.9235]

