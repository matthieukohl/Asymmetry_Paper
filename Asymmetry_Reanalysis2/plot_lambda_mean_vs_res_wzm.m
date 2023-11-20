% plot the mean asymmetry between 925hPA and 200hPA for ncep2, erai, era5 for years 2009-2018 and the
% asymmetry predicted by the moist modal theory vs. moist turbulent theory

close all; clear

% averaging window

latN = 70; latS = 30;
plow = 925; pup = 200;

% ncep2

load('data_ncep2/asymmetry_ncep2_xt_wzm.mat');
load('data_ncep2/asymmetry_ncep2_x_wzm.mat')

level = double(level); lambda_ncep2_x_wzm = double(lambda_ncep2_x_wzm);
lambda_ncep2_xt = double(lambda_ncep2_xt_wzm);

[~,lat30] = min(abs(lat-latS)); [~,latm30] = min(abs(lat+latS));
[~,lat70] = min(abs(lat-latN)); [~,latm70] = min(abs(lat+latN));

[~,p2] = min(abs(level-plow)); [~,p1] = min(abs(level-pup));
dp = level(p1)-level(p2);

lambda_ncep2_NH_x = 1/dp*trapz(level(p2:p1),squeeze(mean(lambda_ncep2_x_wzm(lat70:lat30,p2:p1,:),1)));
lambda_ncep2_SH_x = 1/dp*trapz(level(p2:p1),squeeze(mean(lambda_ncep2_x_wzm(latm30:latm70,p2:p1,:),1)));

lambda_ncep2_NH_xt = 1/dp*trapz(level(p2:p1),squeeze(mean(lambda_ncep2_xt_wzm(lat70:lat30,p2:p1,:),1)));
lambda_ncep2_SH_xt = 1/dp*trapz(level(p2:p1),squeeze(mean(lambda_ncep2_xt_wzm(latm30:latm70,p2:p1,:),1)));

clear('level','lat','p500','lat30','lat70','latm30','latm70',...
    'lambda_ncep2_x','lambda_ncep2_xt');

% erai

load('data_erai/asymmetry_erai_NH_x_wzm.mat'); 
load('data_erai/asymmetry_erai_NH_xt_wzm.mat'); 

level = double(level);

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p2] = min(abs(level-plow)); [~,p1] = min(abs(level-pup));
dp = level(p2)-level(p1);

lambda_erai_NH_x = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_erai_NH_x_wzm(lat70:lat30,p1:p2,:),1)),1);
lambda_erai_NH_xt = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_erai_NH_xt_wzm(lat70:lat30,p1:p2,:),1)));

clear('level','lat','p500','lat30','lat70')

load('data_erai/asymmetry_erai_SH_x_wzm.mat'); 
load('data_erai/asymmetry_erai_SH_xt_wzm.mat');

level = double(level);

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p2] = min(abs(level-plow)); [~,p1] = min(abs(level-pup));
dp = level(p2)-level(p1);

lambda_erai_SH_x = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_erai_SH_x_wzm(latm30:latm70,p1:p2,:),1)));
lambda_erai_SH_xt = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_erai_SH_xt_wzm(latm30:latm70,p1:p2,:),1)));

clear('level','lat','p500','latm30','latm70')

% era5

load('data_era5/asymmetry_era5_NH_x_wzm.mat'); 
load('data_era5/asymmetry_era5_NH_xt_wzm.mat'); 

level = double(level);

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p2] = min(abs(level-plow)); [~,p1] = min(abs(level-pup));
dp = level(p2)-level(p1);

lambda_era5_NH_x = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_era5_NH_x_wzm(lat70:lat30,p1:p2,:),1)));
lambda_era5_NH_xt = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_era5_NH_xt_wzm(lat70:lat30,p1:p2,:),1)));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_wzm.mat'); 
load('data_era5/asymmetry_era5_SH_xt_wzm.mat');

level = double(level);

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p2] = min(abs(level-plow)); [~,p1] = min(abs(level-pup));
dp = level(p2)-level(p1);

lambda_era5_SH_x = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_era5_SH_x_wzm(latm30:latm70,p1:p2,:),1)));
lambda_era5_SH_xt = 1/dp*trapz(level(p1:p2),squeeze(mean(lambda_era5_SH_xt_wzm(latm30:latm70,p1:p2,:),1)));

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','k'); hold on
plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','k','linestyle',':'); hold on;
plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','g'); hold on;
plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
ylabel('\lambda'); ylim([0.5 0.75]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Location','Northwest'); legend boxoff
title(['\lambda (',num2str(plow),'-',num2str(pup),'hPa) wzm: ',num2str(latS),'-',num2str(latN),' NH'])

subplot(1,2,2)

plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','k'); hold on
plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','k','linestyle',':'); hold on;
plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r'); hold on;
plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','g'); hold on;
plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
ylabel('\lambda'); ylim([0.5 0.75])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Location','Northeast'); legend boxoff
title(['\lambda (',num2str(plow),'-',num2str(pup),'hPa) wzm: ',num2str(latS),'-',num2str(latN),' SH'])
