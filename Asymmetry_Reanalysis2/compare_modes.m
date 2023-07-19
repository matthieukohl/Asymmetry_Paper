% compare modal asymmetries

close all; clear;

%load('mode_asymmetry.mat');
load('data_theory/theory_mode_era5_coarse_grained_1p25_res_small_domain.mat')

plot(lambda_NH_modal_theory,'b'); hold on;
plot(lambda_SH_modal_theory,'r');

