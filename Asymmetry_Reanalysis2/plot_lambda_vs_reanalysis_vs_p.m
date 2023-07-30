% plot the lambda from reanalysis vs theory (calculated based on
% time averaged r's and k's

close all; clear;

tic

% averaging window

latN = 70; latS = 30;
lv = 500;

% ncep2

load('data_ncep2/asymmetry_ncep2_xt.mat');
load('data_ncep2/asymmetry_ncep2_x.mat')

[~,lat30] = min(abs(lat-latS)); [~,latm30] = min(abs(lat+latS));
[~,lat70] = min(abs(lat-latN)); [~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

level_ncep2 = level;

lambda_ncep2_NH_x = squeeze(mean(lambda_ncep2_x(lat70:lat30,:,:),1));
lambda_ncep2_SH_x = squeeze(mean(lambda_ncep2_x(latm30:latm70,:,:),1));

lambda_ncep2_NH_xt = squeeze(mean(lambda_ncep2_xt(lat70:lat30,:,:),1));
lambda_ncep2_SH_xt = squeeze(mean(lambda_ncep2_xt(latm30:latm70,:,:),1));

clear('level','lat','p500','lat30','lat70','latm30','latm70',...
    'lambda_ncep2_x','lambda_ncep2_xt');

% erai

load('data_erai/asymmetry_erai_NH_x.mat'); 
load('data_erai/asymmetry_erai_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

level_erai = level;

lambda_erai_NH_x = squeeze(mean(lambda_erai_NH_x(lat70:lat30,:,:),1));
lambda_erai_NH_xt = squeeze(mean(lambda_erai_NH_xt(lat70:lat30,:,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_erai/asymmetry_erai_SH_x.mat'); 
load('data_erai/asymmetry_erai_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

lambda_erai_SH_x = squeeze(mean(lambda_erai_SH_x(latm30:latm70,:,:),1));
lambda_erai_SH_xt = squeeze(mean(lambda_erai_SH_xt(latm30:latm70,:,:),1));

clear('level','lat','p500','latm30','latm70')

% era5

load('data_era5/asymmetry_era5_NH_x.mat'); 
load('data_era5/asymmetry_era5_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
level_era5 = level;

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,:,:),1));
lambda_era5_NH_xt = squeeze(mean(lambda_era5_NH_xt(lat70:lat30,:,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,:,:),1));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,:,:),1));

clear('level','lat','p500','lat30','lat70')

% era5 - coarse

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain.mat');  

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
level_era5 = level;

lambda_era5_NH_x_coarse = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,:,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain.mat'); 

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

lambda_era5_SH_x_coarse = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,:,:),1));

clear('level','lat','p500','lat30','lat70')

% make figure

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)
plot(mean(lambda_ncep2_NH_x(:,[6:8]),2),level_ncep2,'g'); hold on;
plot(mean(lambda_erai_NH_x(:,[6:8]),2),level_erai,'r'); hold on;
plot(mean(lambda_era5_NH_x(:,[6:8]),2),level_era5,'b'); hold on;
plot(mean(lambda_era5_NH_x_coarse(:,[6:8]),2),level_era5,'b--'); hold on;
title('\rm Summer')
set(gca,'Ydir','reverse')
legend('NCEP2','ERAI','ERA5','ERA5 Coarse'); legend boxoff
xlabel('$\lambda$')
ylabel('Pressure (hPa)')

subplot(1,2,2)
plot(mean(lambda_ncep2_NH_x(:,[1,2,12]),2),level_ncep2,'g'); hold on;
plot(mean(lambda_erai_NH_x(:,[1,2,12]),2),level_erai,'r'); hold on;
plot(mean(lambda_era5_NH_x(:,[1,2,12]),2),level_era5,'b'); hold on;
plot(mean(lambda_era5_NH_x_coarse(:,[1,2,12]),2),level_era5,'b--'); hold on;
title('\rm Winter')
set(gca,'Ydir','reverse')
xlabel('$\lambda$')
ylabel('Pressure (hPa)')