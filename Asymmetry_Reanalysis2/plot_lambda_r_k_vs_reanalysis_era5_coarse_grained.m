% compare r and k at 500hpa between ncep2, erai and era5

close all; clear;

% constants for deformation radius

dp = 800*1e2;
f = 1e-4;

% load ncep2 data

%load('data_theory/theory_ncep2_NH_april29th.mat','r_500_NH','k_NH_500','S_500_NH')
%load('data_theory/theory_ncep2_SH_april29th.mat','r_500_SH','k_SH_500','S_500_SH')

load('data_theory/theory_ncep2_NH_july12th_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
load('data_theory/theory_ncep2_SH_july12th_weighted.mat','r_500_SH','k_SH_500','S_500_SH')


r_500_NH_ncep2 = r_500_NH;
r_500_SH_ncep2 = r_500_SH;

% calculate k 

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_erai = LD_NH_500;
LD_SH_500_erai = LD_SH_500;

k_500_NH_ncep2 = k_NH_500.*LD_NH_500_erai;
k_500_SH_ncep2 = k_SH_500.*LD_SH_500_erai;

clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH');

% load data erai

%load('data_theory/theory_erai_NH_april29th.mat','r_500_NH','k_NH_500','S_500_NH')
%load('data_theory/theory_erai_SH_april29th.mat','r_500_SH','k_SH_500','S_500_SH')

load('data_theory/theory_erai_NH_july12th_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
load('data_theory/theory_erai_SH_july12th_weighted.mat','r_500_SH','k_SH_500','S_500_SH')


r_500_NH_erai = r_500_NH;
r_500_SH_erai = r_500_SH;

% calculate k 

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_erai = LD_NH_500;
LD_SH_500_erai = LD_SH_500;

k_500_NH_erai = k_NH_500.*LD_NH_500_erai;
k_500_SH_erai = k_SH_500.*LD_SH_500_erai;

clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH');

% load data era5

% native resolution

% load('data_theory/theory_NH_dec7th.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_SH_dec7th.mat','r_500_SH','k_SH_500','S_500_SH')

% coarse grained

% load('data_theory/theory_era5_NH_coarse_grained_july_5th.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th.mat','r_500_SH','k_SH_500','S_500_SH')

load('data_theory/theory_era5_NH_coarse_grained_res_1p5_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
load('data_theory/theory_era5_SH_coarse_grained_res_1p5_weighted.mat','r_500_SH','k_SH_500','S_500_SH')


% load('data_theory/theory_era5_NH_coarse_grained_july_5th_res_1p75.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th_res_1p75.mat','r_500_SH','k_SH_500','S_500_SH')

r_500_NH_era5 = r_500_NH;
r_500_SH_era5 = r_500_SH;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_era5 = LD_NH_500;
LD_SH_500_era5 = LD_SH_500;

k_500_NH_era5 = k_NH_500.*LD_NH_500_era5;
k_500_SH_era5 = k_SH_500.*LD_SH_500_era5;

r_500_NH_era5_coarse = r_500_NH_era5;
k_500_NH_era5_coarse = k_500_NH_era5;

r_500_SH_era5_coarse = r_500_SH_era5;
k_500_SH_era5_coarse = k_500_SH_era5;

clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH');

% native resolution

% load('data_theory/theory_NH_dec7th.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_SH_dec7th.mat','r_500_SH','k_SH_500','S_500_SH')

load('data_theory/theory_NH_july10th_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
load('data_theory/theory_SH_july10th_weighted.mat','r_500_SH','k_SH_500','S_500_SH')


r_500_NH_era5 = r_500_NH;
r_500_SH_era5 = r_500_SH;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_era5 = LD_NH_500;
LD_SH_500_era5 = LD_SH_500;

k_500_NH_era5 = k_NH_500.*LD_NH_500_era5;
k_500_SH_era5 = k_SH_500.*LD_SH_500_era5;

clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH');


% plot r vs reanalysis

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

figure(1);
x0=10; y0=10; width=1090; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(r_500_NH_ncep2,'g','Linewidth',1.5); hold on;
plot(r_500_NH_erai,'r','Linewidth',1.5); hold on;
plot(r_500_NH_era5_coarse,'m','Linewidth',1.5); hold on;
plot(r_500_NH_era5,'b','Linewidth',1.5); hold on;
legend('NCEP2','ERAI','ERA5 Coarse','ERA5','Location','SouthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Reduction factor r')
ylim([0 0.5])
yticks([0 0.1 0.2 0.3 0.4 0.5])
title('Northern Hemisphere')

subplot(1,2,2)

plot(r_500_SH_ncep2,'g','Linewidth',1.5); hold on;
plot(r_500_SH_erai,'r','Linewidth',1.5); hold on;
plot(r_500_SH_era5_coarse,'b','Linewidth',1.5); hold on;
plot(r_500_SH_era5,'m','Linewidth',1.5); hold on;
legend('NCEP2','ERAI','ERA5 Coarse','ERA5','Location','SouthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Reduction factor r')
ylim([0 0.5])
yticks([0 0.1 0.2 0.3 0.4 0.5])
title('Southern Hemisphere')


% plot k vs reanalysis

figure(2);
x0=10; y0=10; width=1090; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(k_500_NH_era5,'m','Linewidth',1.5); hold on;
plot(k_500_NH_era5_coarse,'b','Linewidth',1.5); hold on;
plot(k_500_NH_erai,'r','Linewidth',1.5); hold on;
plot(k_500_NH_ncep2,'g','Linewidth',1.5); hold on;
legend('ERA5','ERA5 Coarse','ERAI','NCEP2','Location','West'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Wavenumber k')
title('Northern Hemisphere')

subplot(1,2,2)

plot(k_500_SH_era5,'m','Linewidth',1.5); hold on;
plot(k_500_SH_era5_coarse,'b','Linewidth',1.5); hold on;
plot(k_500_SH_erai,'r','Linewidth',1.5); hold on;
plot(k_500_SH_ncep2,'g','Linewidth',1.5); hold on;
legend('ERA5','ERA5 Coarse','ERAI','NCEP2','Location','West'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Wavenumber k')
title('Southern Hemisphere')

% test the theory predictions for reanalysis

% load lambda from reanalysis

% averaging window

latN = 70; latS = 30;
lv = 500;

% ncep2

load('data_ncep2/asymmetry_ncep2_xt.mat');
load('data_ncep2/asymmetry_ncep2_x.mat')

[~,lat30] = min(abs(lat-latS)); [~,latm30] = min(abs(lat+latS));
[~,lat70] = min(abs(lat-latN)); [~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_ncep2_NH_x = squeeze(mean(lambda_ncep2_x(lat70:lat30,p500,:),1));
lambda_ncep2_SH_x = squeeze(mean(lambda_ncep2_x(latm30:latm70,p500,:),1));

lambda_ncep2_NH_xt = squeeze(mean(lambda_ncep2_xt(lat70:lat30,p500,:),1));
lambda_ncep2_SH_xt = squeeze(mean(lambda_ncep2_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70','latm30','latm70',...
    'lambda_ncep2_x','lambda_ncep2_xt');

% erai

load('data_erai/asymmetry_erai_NH_x.mat'); 
load('data_erai/asymmetry_erai_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p500] = min(abs(level-lv));

lambda_erai_NH_x = squeeze(mean(lambda_erai_NH_x(lat70:lat30,p500,:),1));
lambda_erai_NH_xt = squeeze(mean(lambda_erai_NH_xt(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_erai/asymmetry_erai_SH_x.mat'); 
load('data_erai/asymmetry_erai_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_erai_SH_x = squeeze(mean(lambda_erai_SH_x(latm30:latm70,p500,:),1));
lambda_erai_SH_xt = squeeze(mean(lambda_erai_SH_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','latm30','latm70')

% era5

% native grid

load('data_era5/asymmetry_era5_NH_x.mat'); 
load('data_era5/asymmetry_era5_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,p500,:),1));
lambda_era5_NH_xt = squeeze(mean(lambda_era5_NH_xt(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');


[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,p500,:),1));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

% coarse grain

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p25_res.mat'); 
%load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p75_res.mat'); 


[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

lambda_era5_NH_x_coarse_grain = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p25_res.mat');
%load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p75_res.mat');


[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x_coarse_grain = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

% theory

lambda_ncep2_NH_theory = zeros(12,1);
lambda_erai_NH_theory = zeros(12,1);
lambda_era5_NH_theory = zeros(12,1);
lambda_era5_NH_theory_coarse = zeros(12,1);

lambda_ncep2_SH_theory = zeros(12,1);
lambda_erai_SH_theory = zeros(12,1);
lambda_era5_SH_theory = zeros(12,1);
lambda_era5_SH_theory_coarse = zeros(12,1);

for ii = 1:12
    
% % 1d model
% 
lambda_ncep2_NH_theory(ii) = toy_model(r_500_NH_ncep2(ii),k_500_NH_ncep2(ii),1);
lambda_erai_NH_theory(ii) = toy_model(r_500_NH_erai(ii),k_500_NH_erai(ii),1);
lambda_era5_NH_theory(ii) = toy_model(r_500_NH_era5(ii),k_500_NH_era5(ii),1);
lambda_era5_NH_theory_coarse(ii) = toy_model(r_500_NH_era5_coarse(ii),k_500_NH_era5_coarse(ii),1);


lambda_ncep2_SH_theory(ii) = toy_model(r_500_SH_ncep2(ii),k_500_SH_ncep2(ii),1);
lambda_erai_SH_theory(ii) = toy_model(r_500_SH_erai(ii),k_500_SH_erai(ii),1);
lambda_era5_SH_theory(ii) = toy_model(r_500_SH_era5(ii),k_500_SH_era5(ii),1);
lambda_era5_SH_theory_coarse(ii) = toy_model(r_500_SH_era5_coarse(ii),k_500_SH_era5_coarse(ii),1);

% % 2d model
% 
% lambda_ncep2_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_erai(ii),k_500_NH_erai(ii));
% lambda_era5_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_era5(ii),k_500_NH_era5(ii));
% 
% lambda_ncep2_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_erai(ii),k_500_SH_erai(ii));
% lambda_era5_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_era5(ii),k_500_SH_era5(ii));

% use centroid as estimate of ku; use scaling theory to get at kd

% 1d box model

% lambda_ncep2_NH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_NH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_NH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_NH_erai(ii),k_500_NH_erai(ii));
% lambda_era5_NH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_NH_era5(ii),k_500_NH_era5(ii));
% 
% lambda_ncep2_SH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_SH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_SH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_SH_erai(ii),k_500_SH_erai(ii));
% lambda_era5_SH_theory(ii) = toy_model_1d_non_sine_ku_input(r_500_SH_era5(ii),k_500_SH_era5(ii));

% 2d box model 

% lambda_ncep2_NH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_NH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_NH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_NH_erai(ii),k_500_NH_erai(ii));
% lambda_era5_NH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_NH_era5(ii),k_500_NH_era5(ii));
% 
% lambda_ncep2_SH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_SH_ncep2(ii),k_500_NH_ncep2(ii));
% lambda_erai_SH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_SH_erai(ii),k_500_SH_erai(ii));
% lambda_era5_SH_theory(ii) = toy_model_2d_non_sine_ku_input(r_500_SH_era5(ii),k_500_SH_era5(ii));


% 1d and 2d sinusoid models

% x0 = 1;
% ku = k_500_NH_ncep2(ii);
% r = r_500_NH_ncep2(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% K_ncep2(ii) = k;
% 
% lambda_ncep2_NH_theory(ii) = toy_model(r_500_NH_ncep2(ii),k,1);
% 
% x0 = 1;
% ku = k_500_NH_erai(ii);
% r = r_500_NH_erai(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% K_erai(ii) = k;
% 
% lambda_erai_NH_theory(ii) = toy_model(r_500_NH_erai(ii),k,1);
% 
% x0 = 1;
% ku = k_500_NH_era5(ii);
% r = r_500_NH_era5(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% K_era5(ii) = k;
% 
% lambda_era5_NH_theory(ii) = toy_model(r_500_NH_era5(ii),k,1);
% 
% x0 = 1;
% ku = k_500_SH_ncep2(ii);
% r = r_500_SH_ncep2(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% 
% lambda_ncep2_SH_theory(ii) = toy_model(r_500_SH_ncep2(ii),k,1);
% 
% x0 = 1;
% ku = k_500_SH_erai(ii);
% r = r_500_SH_erai(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% 
% lambda_erai_SH_theory(ii) = toy_model(r_500_SH_erai(ii),k,1);
% 
% x0 = 1;
% ku = k_500_SH_era5(ii);
% r = r_500_SH_era5(ii);
% kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);
% k = 2*ku*kd/(ku+kd);
% 
% lambda_era5_SH_theory(ii) = toy_model(r_500_SH_era5(ii),k,1);


% % 2d model
% 
% r = r_500_NH_ncep2(ii); 
% ku = k_500_NH_ncep2(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% K_ncep2(ii) = k;
% 
% 
% lambda_ncep2_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_ncep2(ii),k);
% 
% r = r_500_NH_erai(ii);
% ku = k_500_NH_erai(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% K_erai(ii) = k;
% 
% 
% lambda_erai_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_erai(ii),k);
% 
% r = r_500_NH_era5(ii);
% ku = k_500_NH_era5(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% K_era5(ii) = k;
% 
% 
% lambda_era5_NH_theory(ii) = toy_model_2d_alternative(r_500_NH_era5(ii),k);
% 
% 
% r = r_500_SH_ncep2(ii);
% ku = k_500_SH_ncep2(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% 
% lambda_ncep2_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_ncep2(ii),k);
% 
% 
% r = r_500_SH_erai(ii);
% ku = k_500_SH_erai(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% 
% lambda_erai_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_erai(ii),k);
% 
% 
% r = r_500_SH_era5(ii);
% ku = k_500_SH_era5(ii);
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% 
% lambda_era5_SH_theory(ii) = toy_model_2d_alternative(r_500_SH_era5(ii),k);



end

figure(3);
x0=10; y0=10; width=1090; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(lambda_era5_NH_theory,'m--','Linewidth',1.5); hold on;
plot(lambda_era5_NH_theory_coarse,'b--','Linewidth',1.5); hold on;
plot(lambda_erai_NH_theory,'r--','Linewidth',1.5); hold on;
plot(lambda_ncep2_NH_theory,'g--','Linewidth',1.5); hold on;
plot(lambda_era5_NH_x,'m','Linewidth',1.5); hold on;
plot(lambda_era5_NH_x_coarse_grain,'b','Linewidth',1.5); hold on;
plot(lambda_erai_NH_x,'r','Linewidth',1.5); hold on;
plot(lambda_ncep2_NH_x,'g','Linewidth',1.5); hold on;
legend('ERA5 Theory','ERA5 Theory Coarse','ERAI Theory','NCEP2 Theory','ERA5','ERA5 Coarse','ERAI','NCEP2','Location','NorthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Asymmetry parameter $\lambda$')
title('Northern Hemisphere 1d')
ylim([0.5 0.85])

subplot(1,2,2)
plot(lambda_era5_SH_theory,'m--','Linewidth',1.5); hold on;
plot(lambda_era5_SH_theory_coarse,'b--','Linewidth',1.5); hold on;
plot(lambda_erai_SH_theory,'r--','Linewidth',1.5); hold on;
plot(lambda_ncep2_SH_theory,'g--','Linewidth',1.5); hold on;
plot(lambda_era5_SH_x,'m','Linewidth',1.5); hold on;
plot(lambda_era5_SH_x_coarse_grain,'b','Linewidth',1.5); hold on;
plot(lambda_erai_SH_x,'r','Linewidth',1.5); hold on;
plot(lambda_ncep2_SH_x,'g','Linewidth',1.5); hold on;
legend('ERA5 Theory','ERA5 Theory Coarse','ERAI Theory','NCEP2 Theory','ERA5','ERA5 Coarse','ERAI','NCEP2','Location','NorthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Asymmetry parameter $\lambda$')
title('Southern Hemisphere 1d')
ylim([0.5 0.85])

