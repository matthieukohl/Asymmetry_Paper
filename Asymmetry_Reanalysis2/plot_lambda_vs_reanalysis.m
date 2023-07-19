% plot the lambda from reanalysis vs theory (calculated based on
% time averaged r's and k's

close all; clear;

tic

% averaging window

latN = 70; latS = 30;
lv = 500;

% ncep2

load('data_ncep2/asymmetry_ncep2_xt.mat');
load('data_ncep2/asymmetry_ncep2_x.mat');


[~,lat30] = min(abs(lat-latS)); [~,latm30] = min(abs(lat+latS));
[~,lat70] = min(abs(lat-latN)); [~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

weightN = repmat(cosd(lat(lat70:lat30)),1,1,12);
weightS = repmat(cosd(lat(latm30:latm70)),1,1,12);


lambda_ncep2_NH_x_weighted = squeeze(sum(weightN.*lambda_ncep2_x(lat70:lat30,p500,:),1)./sum(weightN,1));
lambda_ncep2_SH_x_weighted = squeeze(sum(weightS.*lambda_ncep2_x(latm30:latm70,p500,:),1)./sum(weightS,1));

lambda_ncep2_NH_xt_weighted = squeeze(sum(weightN.*lambda_ncep2_xt(lat70:lat30,p500,:),1)./sum(weightN,1));
lambda_ncep2_SH_xt_weighted = squeeze(sum(weightS.*lambda_ncep2_xt(latm30:latm70,p500,:),1)./sum(weightS,1));


lambda_ncep2_NH_x = squeeze(mean(lambda_ncep2_x(lat70:lat30,p500,:),1));
lambda_ncep2_SH_x = squeeze(mean(lambda_ncep2_x(latm30:latm70,p500,:),1));

lambda_ncep2_NH_xt = squeeze(mean(lambda_ncep2_xt(lat70:lat30,p500,:),1));
lambda_ncep2_SH_xt = squeeze(mean(lambda_ncep2_xt(latm30:latm70,p500,:),1));

load('data_ncep2/asymmetry_ncep2_xy_NH.mat')
lambda_ncep2_NH_xy = squeeze(lambda_ncep2_xy(p500,:));

load('data_ncep2/asymmetry_ncep2_xy_SH.mat')
lambda_ncep2_SH_xy = squeeze(lambda_ncep2_xy(p500,:));

clear('level','lat','p500','lat30','lat70','latm30','latm70',...
    'lambda_ncep2_x','lambda_ncep2_xt');

% erai

load('data_erai/asymmetry_erai_NH_x.mat'); 
load('data_erai/asymmetry_erai_NH_xy.mat'); 
load('data_erai/asymmetry_erai_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p500] = min(abs(level-lv));

weightN = repmat(cosd(lat(lat70:lat30)),1,1,12);

lambda_erai_NH_x_weighted = squeeze(sum(weightN.*lambda_erai_NH_x(lat70:lat30,p500,:),1)./sum(weightN,1));

lambda_erai_NH_xt_weighted = squeeze(sum(weightN.*lambda_erai_NH_xt(lat70:lat30,p500,:),1)./sum(weightN,1));


lambda_erai_NH_x = squeeze(mean(lambda_erai_NH_x(lat70:lat30,p500,:),1));
lambda_erai_NH_xy = squeeze(lambda_erai_NH_xy(p500,:));
lambda_erai_NH_xt = squeeze(mean(lambda_erai_NH_xt(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_erai/asymmetry_erai_SH_x.mat'); 
load('data_erai/asymmetry_erai_SH_xy.mat'); 
load('data_erai/asymmetry_erai_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

weightS = repmat(cosd(lat(latm30:latm70)),1,1,12);

lambda_erai_SH_x_weighted = squeeze(sum(weightS.*lambda_erai_SH_x(latm30:latm70,p500,:),1)./sum(weightS,1));
lambda_erai_SH_xt_weighted = squeeze(sum(weightS.*lambda_erai_SH_xt(latm30:latm70,p500,:),1)./sum(weightS,1));


lambda_erai_SH_x = squeeze(mean(lambda_erai_SH_x(latm30:latm70,p500,:),1));
lambda_erai_SH_xy = squeeze(lambda_erai_SH_xy(p500,:));
lambda_erai_SH_xt = squeeze(mean(lambda_erai_SH_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','latm30','latm70')

% era5

load('data_era5/asymmetry_era5_NH_x.mat'); 
load('data_era5/asymmetry_era5_NH_xy_weighted.mat'); 
load('data_era5/asymmetry_era5_NH_xt.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

weightN = repmat(cosd(lat(lat70:lat30)),1,1,12);

lambda_era5_NH_x_weighted = squeeze(sum(weightN.*lambda_era5_NH_x(lat70:lat30,p500,:),1)./sum(weightN,1));

lambda_era5_NH_xt_weighted = squeeze(sum(weightN.*lambda_era5_NH_xt(lat70:lat30,p500,:),1)./sum(weightN,1));


lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,p500,:),1));
lambda_era5_NH_xy = squeeze(lambda_era5_NH_xy(p500,:));
lambda_era5_NH_xt = squeeze(mean(lambda_era5_NH_xt(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xy.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

weightS = repmat(cosd(lat(latm30:latm70)),1,1,12);

lambda_era5_SH_x_weighted = squeeze(sum(weightS.*lambda_era5_SH_x(latm30:latm70,p500,:),1)./sum(weightS,1));
lambda_era5_SH_xt_weighted = squeeze(sum(weightS.*lambda_era5_SH_xt(latm30:latm70,p500,:),1)./sum(weightS,1));


lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,p500,:),1));
lambda_era5_SH_xy = squeeze(lambda_era5_SH_xy(p500,:));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

% % era5 - coarse graining of era5 data
% 
% load('data_era5/asymmetry_era5_NH_x_true_coarse_grain.mat'); 
% load('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p25_res.mat'); 
% 
% 
% [~,lat30] = min(abs(lat-latS)); 
% [~,lat70] = min(abs(lat-latN)); 
% 
% [~,p500] = min(abs(level-lv));
%  
% lambda_era5_NH_x_coarse = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));
% lambda_era5_NH_xy_coarse = squeeze(lambda_era5_NH_xy_coarse_grain(p500,:),1);
% 
% clear('level','lat','p500','lat30','lat70')
% 
% load('data_era5/asymmetry_era5_SH_x_true_coarse_grain.mat'); 
% load('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p25_res.mat'); 
% 
% 
% [~,latm30] = min(abs(lat+latS));
% [~,latm70] = min(abs(lat+latN));
% 
% [~,p500] = min(abs(level-lv));
% 
% lambda_era5_SH_x_coarse = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));
% lambda_era5_SH_xy_coarse = squeeze(lambda_era5_SH_xy_coarse_grain(p500,:));
% 
% 
% clear('level','lat','p500','lat30','lat70')

% era5 - coarse graining of era5 data: 1.5/1.5 resolution

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p5_res.mat'); 
load('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p5_res.mat'); 


[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 

[~,p500] = min(abs(level-lv));

weightN = repmat(cosd(lat(lat70:lat30)),1,1,12);

lambda_era5_NH_x_coarse_1p25_res_weighted = squeeze(sum(weightN.*lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1)./sum(weightN,1));

lambda_era5_NH_x_coarse_1p25_res = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));
lambda_era5_NH_xy_coarse_1p25_res = squeeze(lambda_era5_NH_xy_coarse_grain(p500,:));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p5_res.mat'); 
load('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p5_res.mat'); 


[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

weightS = repmat(cosd(lat(latm30:latm70)),1,1,12);

lambda_era5_SH_x_coarse_1p25_res_weighted = squeeze(sum(weightS.*lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1)./sum(weightS,1));


lambda_era5_SH_x_coarse_1p25_res = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));
lambda_era5_SH_xy_coarse_1p25_res = squeeze(lambda_era5_SH_xy_coarse_grain(p500,:));

clear('level','lat','p500','lat30','lat70')

% make figure

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

h1 = plot(lambda_ncep2_NH_x,'g','Linewidth',1.5); hold on
%plot(lambda_ncep2_NH_xy,'g--','Linewidth',1.5); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%plot(lambda_ncep2_NH_x_weighted,'Linewidth',1.5,'Color','g','linestyle',':'); hold on;

h2 = plot(lambda_erai_NH_x,'r','Linewidth',1.5); hold on;
%plot(lambda_erai_NH_xy,'r--','Linewidth',1.5); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
%plot(lambda_erai_NH_x_weighted,'linewidth',1.5,'Color','r','linestyle',':'); hold on;

h3 = plot(lambda_era5_NH_x,'b--','linewidth',1.5); hold on;
%plot(lambda_era5_NH_xy,'b:','linewidth',1.5); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
%plot(lambda_era5_NH_x_weighted,'linewidth',1.5,'Color','b','Linestyle',':'); hold on;


%h4 = plot(lambda_era5_NH_x_coarse,'b--','linewidth',1.5); hold on;
h5 = plot(lambda_era5_NH_x_coarse_1p25_res,'b','linewidth',1.5); hold on;
%plot(lambda_era5_NH_xy_coarse_1p25_res,'b-.','linewidth',1.5); hold on;
%plot(lambda_era5_NH_x_coarse_1p25_res_weighted,'b-.','linewidth',1.5); hold on;


ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%legend({'NCEP2','ERAI','ERA5','ERA5 Coarse'},'Location','NorthWest'); legend boxoff
l1 = legend([h3,h5,h2,h1],{'ERA5','ERA5 Coarse-Grained','ERAI','NCEP2'},'Location','NorthWest'); legend boxoff
% l1pos = get(l1,'position');
% set(l1,'position',l1pos+[-0.005 0 0 0]);
% l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

h1 = plot(lambda_ncep2_SH_x,'g','Linewidth',1.5); hold on
%plot(lambda_ncep2_SH_xy,'g--','Linewidth',1.5); hold on
%plot(lambda_ncep2_SH_x_weighted,'g--','Linewidth',1.5); hold on


%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
h2 = plot(lambda_erai_SH_x,'r','Linewidth',1.5); hold on;
%plot(lambda_erai_SH_xy,'r--','Linewidth',1.5); hold on;
%plot(lambda_erai_SH_x_weighted,'r--','Linewidth',1.5); hold on;

%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
h3 = plot(lambda_era5_SH_x,'b--','linewidth',1.5); hold on;
%plot(lambda_era5_SH_xy,'b:','linewidth',1.5); hold on;
%plot(lambda_era5_SH_x_weighted,'b:','linewidth',1.5); hold on;

%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
%h4 = plot(lambda_era5_SH_x_coarse,'b--','linewidth',1.5); hold on;
h5 = plot(lambda_era5_SH_x_coarse_1p25_res,'b','linewidth',1.5); hold on;
%plot(lambda_era5_SH_xy_coarse_1p25_res,'b-.','linewidth',1.5); hold on;
%plot(lambda_era5_SH_x_coarse_1p25_res_weighted,'b-.','linewidth',1.5); hold on;


ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1 = legend({'NCEP2','ERAI','ERA5','ERA5 Coarse'},'Location','NorthWest'); legend boxoff
l1 = legend([h3,h5,h2,h1],{'ERA5','ERA5 Coarse-Grained','ERAI','NCEP2'},'Location','NorthWest'); legend boxoff
% l1.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_vs_reanalysis_x_mean','epsc');
