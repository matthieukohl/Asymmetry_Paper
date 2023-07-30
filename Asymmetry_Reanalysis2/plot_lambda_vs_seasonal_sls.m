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


% load r and k values

% load('data_theory/theory_NH_complete2.mat'); 
% load('data_theory/theory_SH_complete2.mat');

load('data_theory/theory_NH_aug2nd.mat'); 
load('data_theory/theory_SH_aug2nd.mat');

% % toy-model
% 
% lambda_NH_theory = zeros(12,1);
% lambda_SH_theory = zeros(12,1);
% 
% for ii = 1:12
%     
% %lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_500(ii)*L_500_NH(ii));
% %lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_500(ii)*L_500_SH(ii));
% lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_500(ii)*1e6/sqrt(2));
% lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_500(ii)*1e6/sqrt(2));
% end

% Paul's method - nondimensionalize using real N2
load('data_theory/theory_NH_oct13th.mat'); 
load('data_theory/theory_SH_oct13th.mat');

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

% toy-model

lambda_NH_theory = zeros(12,1);
lambda_SH_theory = zeros(12,1);

for ii = 1:12
    
%lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_500(ii)*L_500_NH(ii));
%lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_500(ii)*L_500_SH(ii));
lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
%lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_500(ii)*0.5*1e6/sqrt(2));
%lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_500(ii)*0.5*1e6/sqrt(2));
end



% modal theory

lambda_NH_modal_theory = zeros(12,1);
lambda_SH_modal_theory = zeros(12,1);

% for ii = 1:12
%     
% lambda_NH_modal_theory(ii) = mode_asymmetry(r_500_NH(ii));
% lambda_SH_modal_theory(ii) = mode_asymmetry(r_500_SH(ii));
% end

%save('mode_asymmetry_bigger_domain.mat','lambda_NH_modal_theory','lambda_SH_modal_theory');
load('mode_asymmetry_bigger_domain.mat');
load('mode_asymmetry.mat','lambda_SH_modal_theory');


% figure for SLS:: just ERA5
plot(lambda_era5_NH_x,'b','Linewidth',1.5); hold on;

ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylim([0.5 0.75])
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/lambda_seasonal_cycle','epsc');


% all reanalyses 

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
%plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
%plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1pos = get(l1,'position');
%set(l1,'position',l1pos+[-0.005 0 0 0]);
%l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
%title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')
ylim([0.5 0.75])

subplot(1,2,2)

plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r','LineStyle','--'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
%plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
%plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','Location','Northwest'); legend boxoff
%l1.FontSize = 9;
%title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')
ylim([0.5 0.75])


saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/lambda_seasonal_cycle_all_reanalyses','epsc');

hfig = figure(2);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

%plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1pos = get(l1,'position');
%set(l1,'position',l1pos+[-0.005 0 0 0]);
%l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
%title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')
ylim([0.5 0.85])

subplot(1,2,2)

%plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r','LineStyle','--'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','Location','Northwest'); legend boxoff
%l1.FontSize = 9;
%title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')
ylim([0.5 0.85])


saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/lambda_seasonal_cycle_vs_theory','epsc');





% make figure

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','b','LineStyle',':'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','b','LineStyle','--'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_NH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[-0.005 0 0 0]);
l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','b','LineStyle',':'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','b','LineStyle','--'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_SH_theory,'linewidth',1.5,'Color','k'); hold on;
plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
l1.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_reanalysis_vs_theory_new_all_reanalyses','epsc');

%saveas(gcf,[pwd,'/AOFD2022/lambda_season_vs_theory_all_reanalysis'],'epsc');
%saveas(gcf,[pwd,'/AOFD2022/lambda_season_vs_theory_only_era5'],'epsc');

figure 
x0=10; y0=10; width=8/7*450; height=8/7*300; set(gcf,'position',[x0,y0,width,height])
plot(r_500_NH,'k','Linewidth',1.5); hold on;
plot(r_500_SH,'k--','Linewidth',1.5)
legend('NH','SH','Location','SouthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Reduction factor r')
yticks([0.1 0.2 0.3 0.4 0.5])
title('Seasonal Cycle of the Reduction Factor r (500hPa)')
%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/r_seasonal_cycle','epsc');

figure 
x0=10; y0=10; width=8/7*450; height=8/7*300; set(gcf,'position',[x0,y0,width,height])
plot(k_NH_500.*LD_NH_500,'k','Linewidth',1.5); hold on;
plot(k_SH_500.*LD_SH_500,'k--','Linewidth',1.5)
legend('NH','SH','Location','NorthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Reduction factor r')
%yticks([0.1 0.2 0.3 0.4 0.5])
title('Seasonal Cycle of the Wavenumber k (500hPa)')