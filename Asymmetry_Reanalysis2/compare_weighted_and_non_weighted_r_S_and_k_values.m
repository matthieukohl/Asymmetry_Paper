% compare cos-weighted r, S and k values

close all; clear;

% non weighted

% % era5
% load('data_theory/theory_NH_dec7th.mat');
% load('data_theory/theory_SH_dec7th.mat');

load('data_theory/theory_era5_NH_coarse_grained_july_5th.mat');
load('data_theory/theory_era5_SH_coarse_grained_july_5th.mat');

% % erai
% load('data_theory/theory_erai_NH_april29th.mat');
% load('data_theory/theory_erai_SH_april29th.mat');

% % ncep2
% load('data_theory/theory_ncep2_NH_april29th.mat');
% load('data_theory/theory_ncep2_SH_april29th.mat');

r_500_NH_nw = r_500_NH;
r_500_SH_nw = r_500_SH;

S_500_NH_nw = S_500_NH;
S_500_SH_nw = S_500_SH;

k_500_NH_nw = k_NH_500;
k_500_SH_nw = k_SH_500;

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500_nw = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500_nw = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

% weighted

% era5
% load('data_theory/theory_NH_july10th_weighted.mat');
% load('data_theory/theory_SH_july10th_weighted.mat');

load('data_theory/theory_era5_NH_coarse_grained_july_12th_res_1p25_weighted.mat');
load('data_theory/theory_era5_SH_coarse_grained_july_12th_res_1p25_weighted.mat');



% erai

% load('data_theory/theory_erai_NH_july12th_weighted.mat');
% load('data_theory/theory_erai_SH_july12th_weighted.mat');


%ncep2
%load('data_theory/theory_ncep2_NH_july12th_weighted.mat');
%load('data_theory/theory_ncep2_SH_july12th_weighted.mat');


r_500_NH_w = r_500_NH;
r_500_SH_w = r_500_SH;

S_500_NH_w = S_500_NH;
S_500_SH_w = S_500_SH;

k_500_NH_w = k_NH_500;
k_500_SH_w = k_SH_500;

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500_w = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500_w = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height



% make plots

figure

plot(r_500_NH_nw,'b'); hold on;
plot(r_500_NH_w,'b--'); hold on;
plot(r_500_SH_nw,'r'); hold on;
plot(r_500_SH_w,'r--');
ylabel('Reduction factor r')
set(gca,'FontSize',12);
legend('NH','NH weighted','SH','SH weighted','Location','SouthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)

figure

plot(S_500_NH_nw,'b'); hold on;
plot(S_500_NH_w,'b--'); hold on;
plot(S_500_SH_nw,'r'); hold on;
plot(S_500_SH_w,'r--');
set(gca,'FontSize',12);
ylabel('Static Stability')
legend('NH','NH weighted','SH','SH weighted'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)

figure

plot(k_500_NH_nw,'b'); hold on;
plot(k_500_NH_w,'b--'); hold on;
plot(k_500_SH_nw,'r'); hold on;
plot(k_500_SH_w,'r--');
set(gca,'FontSize',12);
ylabel('Wavenumber k')
legend('NH','NH weighted','SH','SH weighted'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)

figure

plot(k_500_NH_nw.*LD_NH_500_nw,'b'); hold on;
plot(k_500_NH_w.*LD_NH_500_w,'b--'); hold on;
plot(k_500_SH_nw.*LD_SH_500_nw,'r'); hold on;
plot(k_500_SH_w.*LD_SH_500_nw,'r--');
set(gca,'FontSize',12);
ylabel('Wavenumber k (nondimensionalized)')
legend('NH','NH weighted','SH','SH weighted'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)




