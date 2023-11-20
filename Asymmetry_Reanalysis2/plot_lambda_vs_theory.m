% plot the asymmetry at 500hPA for ncep2, erai, era5 for years 2009-2018 and the
% asymmetry predicted by the moist modal theory vs. moist turbulent theory

close all; clear

% averaging window

latN = 60; latS = 40;
lv =500;

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

% load theory

% load('data_theory/theory_NH.mat'); 
% load('data_theory/theory_SH.mat');
load('data_theory/theory_NH_complete2.mat'); 
load('data_theory/theory_SH_complete2.mat');
% load('data_theory/theory_NH_randn.mat'); 
% load('data_theory/theory_SH_randn.mat');
% load('data_theory/theory_NH_sin_k_rhs.mat');
% load('data_theory/theory_SH_sin_k_rhs.mat');
% % 
% 
% lambda_NH_theory_mid = lambda_NH_theory_500;
% lambda_SH_theory_mid = lambda_SH_theory_500;
% r_mid_NH = r_500_NH;
% r_mid_SH = r_500_SH;

lambda_NH_theory = lambda_NH_theory_mid;
lambda_SH_theory = lambda_SH_theory_mid;

load('data_theory/modal_theory_nh_N800_L32pi.mat');
lambda_NH_modal_theory = lambda; clear('lambda');
lambda_NH_modal_ZG = 1./(1+sqrt(r_mid_NH));
lambda_NH_modal_PG = 1./(1+r_mid_NH);
load('data_theory/modal_theory_sh_N800_L32pi.mat');
lambda_SH_modal_theory = lambda; clear('lambda');
lambda_SH_modal_ZG = 1./(1+sqrt(r_mid_SH));
lambda_SH_modal_PG = 1./(1+r_mid_SH);


% hfig = figure(1);
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 2 1])
% 
% subplot(1,2,1)
% 
% plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','m'); hold on
% plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
% plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
% plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
% plot(lambda_era5_NH_x,'linewidth',1.5,'Color','g'); hold on;
% plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
% plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
% plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle',':');hold on
% plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
% plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
% ylabel('\lambda'); ylim([0.5 1.0]) 
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',10)
% legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northwest'); legend boxoff
% %legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
% title(['\lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' NH'])
% set(gca,'box','off')
% 
% subplot(1,2,2)
% 
% plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','m'); hold on
% plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
% plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r'); hold on;
% plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
% plot(lambda_era5_SH_x,'linewidth',1.5,'Color','g'); hold on;
% plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
% plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
% plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle',':'); hold on;
% plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
% plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
% ylabel('\lambda'); ylim([0.5 1.0])
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',10)
% %legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
% title(['\lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' SH'])
% set(gca,'box','off')
% saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_vs_seasonal_randn','epsc');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

%plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','m'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter \lambda'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['\rm \lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' NH'])
set(gca,'box','off')

subplot(1,2,2)

%plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','m'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter \lambda'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
title(['\rm \lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' SH'])
set(gca,'box','off')
saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_reanalysis_vs_theory','epsc');
%saveas(gcf,[pwd,'/AOFD2022/lambda_season_vs_theory_all_reanalysis'],'epsc');
%saveas(gcf,[pwd,'/AOFD2022/lambda_season_vs_theory_only_era5'],'epsc');




hfig = figure(2);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','m'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','g'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle',':');hold on
plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('\lambda'); ylim([0.5 1.0]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('NCEP2(x)','ERAI(x)','ERA5(x)','Turb','Modal','ZG','PG','Location','Northwest'); legend boxoff
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['\lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' NH'])
set(gca,'box','off')

subplot(1,2,2)

plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','m'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','g'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle',':'); hold on;
plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('\lambda'); ylim([0.5 1.0])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
title(['\lambda (',num2str(lv),'hPa): ',num2str(latS),'-',num2str(latN),' SH'])
set(gca,'box','off')
saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_reanalysis_vs_theory','epsc');



hfig = figure(2);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)
plot(r_mid_NH,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('r')
title('Northern Hemisphere')

subplot(1,2,2)
plot(r_mid_SH,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('r')
title('Southern Hemisphere')

hfig = figure(3);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)
plot(L_mid_NH/1000,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('L_D(km)')
title('Northern Hemisphere')

subplot(1,2,2)
plot(-L_mid_SH/1000,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('L_D(km)')
title('Southern Hemisphere')

hfig = figure(4);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)
plot(k_NH_mid.*L_mid_NH,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('k(ND)')
title('Northern Hemisphere')

subplot(1,2,2)
plot(-k_SH_mid.*L_mid_SH,'linewidth',1.5)
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
ylabel('k(ND)')
title('Southern Hemisphere')