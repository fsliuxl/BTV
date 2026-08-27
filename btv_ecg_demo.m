clc; clear;
close all;
seed = 2025;
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
addpath(genpath('MyCode'));
addpath(genpath('quality asses'));
%% BTV for signal recovery
% Parameters
var_struct = load(strcat(pwd,'\fecg.mat'));
name_cell = fieldnames(var_struct);
x_ori = double(getfield(var_struct,char(name_cell)));
x_ori = x_ori(1:512,2);
x_ori = (x_ori - mean(x_ori))/max(abs(x_ori));
%% Sampling
n      = length(x_ori);
sr     = .5;
m      = ceil(sr * n); % sampling number
Phi    = randn(m,n); % measurement matrix
b      = Phi*x_ori; % observed signal
%% BTV for signal recovery
lambda      = 1; % regularization parameter
block_size  = 7; % block size
% the above two parameters need to be well tuned in experiments.
opts.lambda = lambda;
opts.block_size = block_size;
fecg_btv = btv_admm(b,Phi,opts);
[myRelErr(real(x_ori),fecg_btv),mysnr(real(x_ori),fecg_btv)]

% The result is [0.0714, 22.9235]


