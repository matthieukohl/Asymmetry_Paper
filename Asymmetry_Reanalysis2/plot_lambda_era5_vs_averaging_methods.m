% plot lambda era5 vs. different averaging methods

close all; clear;

% averaging window

latN = 70; latS = 30;
lv = 500;

% era5

load('data_era5/asymmetry_era5_NH_x.mat'); 
load('data_era5/asymmetry_era5_NH_xy.mat'); 
load('data_era5/asymmetry_era5_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,p500,:),1));
lambda_era5_NH_xy = squeeze(lambda_era5_NH_xy(p500,:));
lambda_era5_NH_xt = squeeze(mean(lambda_era5_NH_xt(lat70:lat30,p500,:),1));

% without zonal mean (wzm)

load('data_era5/asymmetry_era5_NH_x_wzm.mat'); 
load('data_era5/asymmetry_era5_NH_xy_wzm.mat'); 
load('data_era5/asymmetry_era5_NH_xt_wzm.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

lambda_era5_NH_x_wzm = squeeze(mean(lambda_era5_NH_x_wzm(lat70:lat30,p500,:),1));
lambda_era5_NH_xy_wzm = squeeze(lambda_era5_NH_xy_wzm(p500,:));
lambda_era5_NH_xt_wzm = squeeze(mean(lambda_era5_NH_xt_wzm(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xy.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,p500,:),1));
lambda_era5_SH_xy = squeeze(lambda_era5_SH_xy(p500,:));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,p500,:),1));

load('data_era5/asymmetry_era5_SH_x_wzm.mat'); 
load('data_era5/asymmetry_era5_SH_xy_wzm.mat'); 
load('data_era5/asymmetry_era5_SH_xt_wzm.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x_wzm = squeeze(mean(lambda_era5_SH_x_wzm(latm30:latm70,p500,:),1));
lambda_era5_SH_xy_wzm = squeeze(lambda_era5_SH_xy_wzm(p500,:));
lambda_era5_SH_xt_wzm = squeeze(mean(lambda_era5_SH_xt_wzm(latm30:latm70,p500,:),1));


clear('level','lat','p500','lat30','lat70')

% era5 - coarse graining of era5 data

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p500] = min(abs(level-lv));
 
lambda_era5_NH_x_coarse = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain.mat'); 

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x_coarse = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

% era5 - coarse graining of era5 data: 1.25/1.25 resolution

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p25_res.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p500] = min(abs(level-lv));
 
lambda_era5_NH_x_coarse_1p25_res = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p25_res.mat'); 

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x_coarse_1p25_res = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

% make figure

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

h1 = plot(lambda_era5_NH_x,'g','Linewidth',1.5); hold on
h2 = plot(lambda_era5_NH_x_wzm,'g--','Linewidth',1.5); hold on
h3 = plot(lambda_era5_NH_xy,'r','Linewidth',1.5); hold on;
h4 = plot(lambda_era5_NH_xy_wzm,'r--','Linewidth',1.5); hold on;
h5 = plot(lambda_era5_NH_xt,'b','linewidth',1.5); hold on;
h6 = plot(lambda_era5_NH_xt_wzm,'b--','linewidth',1.5); hold on;

ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('x-average','xy-average','xt-average','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%legend({'NCEP2','ERAI','ERA5','ERA5 Coarse'},'Location','NorthWest'); legend boxoff
l1 = legend([h1,h2,h3,h4,h5,h6],{'x-average','x-average wzm','xy-average','xy-average wzm','xt-average','xt-average wzm'},'Location','NorthWest'); legend boxoff
% l1pos = get(l1,'position');
% set(l1,'position',l1pos+[-0.005 0 0 0]);
% l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

h1 = plot(lambda_era5_SH_x,'g','Linewidth',1.5); hold on
h2 = plot(lambda_era5_SH_x_wzm,'g--','Linewidth',1.5); hold on

h3 = plot(lambda_era5_SH_xy,'r','Linewidth',1.5); hold on;
h4 = plot(lambda_era5_SH_xy_wzm,'r--','Linewidth',1.5); hold on;

h5 = plot(lambda_era5_SH_xt,'b','linewidth',1.5); hold on;
h6 = plot(lambda_era5_SH_xt_wzm,'b--','linewidth',1.5); hold on;


ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('x-average','xy-average','xt-average'); legend boxoff
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1 = legend({'NCEP2','ERAI','ERA5','ERA5 Coarse'},'Location','NorthWest'); legend boxoff
l1 = legend([h1,h2,h3,h4,h5,h6],{'x-average','x-average wzm','xy-average','xy-average wzm','xt-average','xt-average wzm'},'Location','NorthWest'); legend boxoff
% l1.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_era5_vs_averaging_methods','epsc');
