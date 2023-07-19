% main mode calculation

% load data era5 coarse grained

% load('data_theory/theory_era5_NH_coarse_grained_july_5th.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th.mat','r_500_SH','k_SH_500','S_500_SH')

%load('data_theory/theory_NH_dec7th.mat')
%load('data_theory/theory_SH_dec7th.mat')

% load('data_theory/theory_era5_NH_coarse_grained_july_12th_res_1p25_weighted.mat');
% load('data_theory/theory_era5_SH_coarse_grained_july_12th_res_1p25_weighted.mat');

load('data_theory/theory_era5_NH_coarse_grained_res_1p5_weighted.mat');
load('data_theory/theory_era5_SH_coarse_grained_res_1p5_weighted.mat');


L = 32*pi;
N = 800;

lambda_NH_modal_theory = zeros(12,1);
lambda_SH_modal_theory = zeros(12,1);
w_final_NH = zeros(N+1,12);
w_final_SH = zeros(N+1,12);

for ii = 1:12
ii
    
[lambda_NH_modal_theory(ii),w_final_NH(:,ii)] = mode_asymmetry(r_500_NH(ii),N,L);
[lambda_SH_modal_theory(ii),w_final_SH(:,ii)] = mode_asymmetry(r_500_SH(ii),N,L);
end


%save('data_theory/theory_mode_era5_coarse_grained_res_1p75_weighted.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory');
%save('data_theory/theory_mode_era5_coarse_grained_res_1p75_weighted_bigger_domain.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory');

save('data_theory/theory_mode_era5_coarse_grained_res_1p5_weighted_even_bigger_domain.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory','w_final_NH','w_final_SH','N','L');

%save('data_theory/theory_mode_era5_coarse_grained_1p25_res_small_domain.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory');

% load

% close all; clear
% 
% load('data_theory/theory_mode_era5_coarse_grained_res_1p5_weighted.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory','w_final_NH','w_final_SH','N','L');
% lambda_NH = lambda_NH_modal_theory;
% lambda_SH = lambda_SH_modal_theory;
% w_NH = w_final_NH;
% w_SH = w_final_SH;
% load('data_theory/theory_mode_era5_coarse_grained_res_1p5_weighted_bigger_domain.mat','r_500_NH','lambda_NH_modal_theory','lambda_SH_modal_theory','w_final_NH','w_final_SH','N','L');
% lambda_NH_bigger_domain = lambda_NH_modal_theory;
% lambda_SH_bigger_domain = lambda_SH_modal_theory;
% w_NH_bigger_domain = w_final_NH;
% w_SH_bigger_domain = w_final_SH;




