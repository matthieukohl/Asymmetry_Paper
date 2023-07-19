% plot the asymmetry at 500hPA for ncep2, erai, era5 for years 2009-2018 and the
% asymmetry predicted by the moist modal theory vs. moist turbulent theory

close all; clear

% averaging window

latN = 70; latS = 30;

% ncep2

load('data_ncep2/asymmetry_ncep2_xt.mat');
load('data_ncep2/asymmetry_ncep2_x.mat')

level_ncep2 = level;

[~,lat30] = min(abs(lat-latS)); [~,latm30] = min(abs(lat+latS));
[~,lat70] = min(abs(lat-latN)); [~,latm70] = min(abs(lat+latN));

lambda_ncep2_NH_x = squeeze(mean(lambda_ncep2_x(lat70:lat30,:,:),1));
lambda_ncep2_SH_x = squeeze(mean(lambda_ncep2_x(latm30:latm70,:,:),1));

lambda_ncep2_NH_xt = squeeze(mean(lambda_ncep2_xt(lat70:lat30,:,:),1));
lambda_ncep2_SH_xt = squeeze(mean(lambda_ncep2_xt(latm30:latm70,:,:),1));

clear('level','lat','p500','lat30','lat70','latm30','latm70',...
    'lambda_ncep2_x','lambda_ncep2_xt');

% erai

load('data_erai/asymmetry_erai_NH_x.mat'); 
load('data_erai/asymmetry_erai_NH_xt.mat'); 

level_erai = level;

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

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

level_era5 = level;

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,:,:),1));
lambda_era5_NH_xt = squeeze(mean(lambda_era5_NH_xt(lat70:lat30,:,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,:,:),1));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,:,:),1));

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(mean(lambda_ncep2_NH_x(2:end,6:8),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m'); hold on
plot(mean(lambda_ncep2_NH_x(2:end,[12,1,2]),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(mean(lambda_erai_NH_x(1:end-3,6:8),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r'); hold on
plot(mean(lambda_erai_NH_x(1:end-3,[12,1,2]),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(mean(lambda_era5_NH_x(1:end-3,6:8),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g'); hold on
plot(mean(lambda_era5_NH_x(1:end-3,[12,1,2]),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g','linestyle',':'); hold on;
ylabel('p(hPA)'); ylim([0,925]); set(gca,'Ydir','reverse');
xlabel('\lambda');
set(gca,'Fontsize',10)
legend('NCEP2(JJA)','NCEP2(DJF)','ERAI(JJA)','ERAI(DJF)','ERA5(JJA)','ERA5(DJF)','Location','West'); legend boxoff
title(['\lambda(p) x-average: ',num2str(latS),'-',num2str(latN),' NH'])

subplot(1,2,2)

plot(mean(lambda_ncep2_SH_x(2:end,6:8),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m'); hold on
plot(mean(lambda_ncep2_SH_x(2:end,[12,1,2]),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(mean(lambda_erai_SH_x(1:end-3,6:8),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r'); hold on
plot(mean(lambda_erai_SH_x(1:end-3,[12,1,2]),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(mean(lambda_era5_SH_x(1:end-3,6:8),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g'); hold on
plot(mean(lambda_era5_SH_x(1:end-3,[12,1,2]),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g','linestyle',':'); hold on;
ylabel('p(hPA)'); ylim([0,925]); set(gca,'Ydir','reverse');
xlabel('\lambda');
set(gca,'Fontsize',10)
legend('NCEP2(JJA)','NCEP2(DJF)','ERAI(JJA)','ERAI(DJF)','ERA5(JJA)','ERA5(DJF)','Location','West'); legend boxoff
title(['\lambda(p) x-average: ',num2str(latS),'-',num2str(latN),' SH'])

hfig = figure(2);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(mean(lambda_ncep2_NH_xt(2:end,6:8),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m'); hold on
plot(mean(lambda_ncep2_NH_xt(2:end,[12,1,2]),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(mean(lambda_erai_NH_xt(1:end-3,6:8),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r'); hold on
plot(mean(lambda_erai_NH_xt(1:end-3,[12,1,2]),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(mean(lambda_era5_NH_xt(1:end-3,6:8),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g'); hold on
plot(mean(lambda_era5_NH_xt(1:end-3,[12,1,2]),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g','linestyle',':'); hold on;
ylabel('p(hPA)'); ylim([0,925]); set(gca,'Ydir','reverse');
xlabel('\lambda');
set(gca,'Fontsize',10)
legend('NCEP2(JJA)','NCEP2(DJF)','ERAI(JJA)','ERAI(DJF)','ERA5(JJA)','ERA5(DJF)','Location','West'); legend boxoff
title(['\lambda(p) xt-average: ',num2str(latS),'-',num2str(latN),' NH'])

subplot(1,2,2)

plot(mean(lambda_ncep2_SH_xt(2:end,6:8),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m'); hold on
plot(mean(lambda_ncep2_SH_xt(2:end,[12,1,2]),2),level_ncep2(2:end),'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(mean(lambda_erai_SH_xt(1:end-3,6:8),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r'); hold on
plot(mean(lambda_erai_SH_xt(1:end-3,[12,1,2]),2),level_erai(1:end-3),'Linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(mean(lambda_era5_SH_xt(1:end-3,6:8),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g'); hold on
plot(mean(lambda_era5_SH_xt(1:end-3,[12,1,2]),2),level_era5(1:end-3),'Linewidth',1.5,'Color','g','linestyle',':'); hold on;
ylabel('p(hPA)'); ylim([0,925]); set(gca,'Ydir','reverse');
xlabel('\lambda');
set(gca,'Fontsize',10)
legend('NCEP2(JJA)','NCEP2(DJF)','ERAI(JJA)','ERAI(DJF)','ERA5(JJA)','ERA5(DJF)','Location','West'); legend boxoff
title(['\lambda(p) xt-average: ',num2str(latS),'-',num2str(latN),' SH'])


