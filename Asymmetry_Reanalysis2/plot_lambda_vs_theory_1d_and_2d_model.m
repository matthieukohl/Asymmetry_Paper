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
% load('data_theory/theory_NH_oct13th.mat'); 
% load('data_theory/theory_SH_oct13th.mat');

load('data_theory/theory_NH_dec7th.mat'); 
load('data_theory/theory_SH_dec7th.mat');

dp = 800*1e2;
f = 1e-4;

% LD_NH_500  = 1e6/(2*sqrt(2))*ones(12,1);
% LD_SH_500 = LD_NH_500;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

% load the deformation radius from N^2/f^2 w_pp / w calculations
% load('data_theory/deformation_radius_era5_NH_april_18th.mat');
% load('data_theory/deformation_radius_era5_SH_april_18th.mat');

% load('data_theory/deformation_radius_era5_NH_april_18th_500hPa.mat');
% load('data_theory/deformation_radius_era5_SH_april_18th_500hPa.mat');

[a,p500] = min(abs(level-500*1e2));

% LD_NH_500 = squeeze(Ld_NH(p500,:));
% LD_SH_500 = squeeze(Ld_SH(p500,:));

% toy-model

lambda_NH_theory_1d = zeros(12,1);
lambda_NH_theory_2d = zeros(12,1);
lambda_NH_theory_2d_refined = zeros(12,1);
lambda_SH_theory_1d = zeros(12,1);
lambda_SH_theory_2d = zeros(12,1);
lambda_SH_theory_2d_refined = zeros(12,10);


for ii = 1:12
    
%lambda_NH_theory_1d(ii) = toy_model_1d_non_sine(r_500_NH(ii),1.71);
%lambda_NH_theory_2d(ii) = toy_model_2d_non_sine(r_500_NH(ii),1.71);
%lambda_NH_theory_2d(ii) = toy_model_2d_non_sine_refined(r_500_NH(ii),1.71);


%lambda_SH_theory_1d(ii) = toy_model_1d_non_sine(r_500_SH(ii),1.71);
%lambda_SH_theory_2d(ii) = toy_model_2d_non_sine(r_500_SH(ii),1.71);
%lambda_SH_theory_2d(ii) = toy_model_2d_non_sine_refined(r_500_SH(ii),1.71);


% use kd from w-spectrum

%lambda_NH_theory_1d(ii) = toy_model_1d_non_sine(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
%lambda_NH_theory_2d(ii) = toy_model_2d_non_sine(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
%lambda_NH_theory_2d(ii) = toy_model_2d_non_sine_refined(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
%lambda_NH_theory_2d(ii) = toy_model_2d_alternative(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
% %lambda_NH_theory_2d_refined(ii) = toy_model_2d(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
% 
% 
%lambda_SH_theory_1d(ii) = toy_model_1d_non_sine(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
%lambda_SH_theory_2d(ii) = toy_model_2d_non_sine(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
%lambda_SH_theory_2d(ii) = toy_model_2d_non_sine_refined(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
%lambda_SH_theory_2d(ii) = toy_model_2d_alternative(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
% %lambda_SH_theory_2d_refined(ii) = toy_model_2d(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));

% assume that centroid of rhs/(1+k^2)^2 spectrum yields kd rather than k 

lambda_NH_theory_1d(ii) = toy_model_1d_non_sine(r_500_NH(ii),0.77);
%lambda_NH_theory_2d(ii) = toy_model_2d_non_sine(r_500_NH(ii),0.77);

r = r_500_NH(ii);

kd = 0.77;

kd2 = kd^2;

ku2 = (-1+sqrt(1+8*r.*kd2.*(1+2*kd2)))./(4*r); % alternative definition

ku = sqrt(ku2);

k = 2*ku*kd/(ku+kd);

lambda_NH_theory_2d(ii) = toy_model_2d_alternative(r_500_NH(ii),k);


lambda_SH_theory_1d(ii) = toy_model_1d_non_sine(r_500_SH(ii),0.77);
%lambda_SH_theory_2d(ii) = toy_model_2d_non_sine(r_500_SH(ii),0.77);

r = r_500_SH(ii);

kd = 0.77;

kd2 = kd^2;

ku2 = (-1+sqrt(1+8*r.*kd2.*(1+2*kd2)))./(4*r); % alternative definition

ku = sqrt(ku2);

k = 2*ku*kd/(ku+kd);

lambda_SH_theory_2d(ii) = toy_model_2d_alternative(r_500_SH(ii),k);



% assume that centroid of w spectrum yields ku rather than k 

% lambda_NH_theory_1d(ii) = toy_model_1d_non_sine_ku_input(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));
% lambda_NH_theory_2d(ii) = toy_model_2d_non_sine_ku_input(r_500_NH(ii),k_NH_500(ii)*LD_NH_500(ii));

% calculate k for sinusoid model from 2/k = 1/ku + 1/kd 
% where ku is inputted and kd is calculated from the scaling theory

% r = r_500_NH(ii);
% 
% ku = k_NH_500(ii)*LD_NH_500(ii);
% 
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% lambda_NH_theory_2d(ii) = toy_model_2d_alternative(r_500_NH(ii),k);


% lambda_SH_theory_1d(ii) = toy_model_1d_non_sine_ku_input(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));
% lambda_SH_theory_2d(ii) = toy_model_2d_non_sine_ku_input(r_500_SH(ii),k_SH_500(ii)*LD_SH_500(ii));

% r = r_500_SH(ii);
% 
% ku = k_SH_500(ii)*LD_SH_500(ii);
% 
% ku2 = ku^2;
% kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));
% 
% kd = sqrt(kd2);
% 
% k = 2*ku*kd/(ku+kd);
% 
% lambda_SH_theory_2d(ii) = toy_model_2d_alternative(r_500_SH(ii),k);
% 




end


% modal theory

lambda_NH_modal_theory = zeros(12,1);
lambda_SH_modal_theory = zeros(12,1);

% need to look into this again; for asymmetry.mat lambda mode is not as
% smooth over the seasonal cycle as r itself. seems that there is a
% jump of lambda for certain r values when an additional 
% updraft/downdraft fits into the box; running on bigger domain helps
% but then the SH, for which r is larger, looks non-smooth

load('mode_asymmetry_bigger_domain.mat');
load('mode_asymmetry.mat','lambda_SH_modal_theory');
%load('mode_asymmetry.mat');


% for ii = 1:12
%     
% lambda_NH_modal_theory(ii) = mode_asymmetry(r_500_NH(ii));
% lambda_SH_modal_theory(ii) = mode_asymmetry(r_500_SH(ii));
% end

% % figure for SLS
% plot(lambda_era5_NH_x,'b','Linewidth',1.5)
% ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% set(gca,'box','off')
% ylim([0.5 0.75])
% saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/lambda_seasonal_cycle','epsc');

% make figure

% set interpreter to latex
set(groot,'DefaultTextInterpreter','latex')
set(groot,'DefaultLegendInterpreter','latex')

hfig = figure(1);
%pos = get(hfig,'position');
%set(hfig,'position',pos.*[.5 1 2 1])
x0=10; y0=10; width=1300; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

%h1 = plot(lambda_ncep2_NH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_NH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%h2 = plot(lambda_erai_NH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_NH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
h3 = plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_NH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
h5 = plot(lambda_NH_theory_1d,'linewidth',1.5,'Color','k'); hold on;
h6 = plot(lambda_NH_theory_2d,'k:','linewidth',1.5); hold on;
%h6 = plot(lambda_NH_theory_2d_refined,'linewidth',1.5,'Color','k','Linestyle',':');hold on
h8 = plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_NH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
%l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1 = legend([h7,h6,h5,h3,h4,h2,h1],{'1-D Modal Theory','1-D Full Spectrum','1-D Toy Model','ERA5','ERA5 Corase','ERAI','NCEP2'},'Location','NorthWest'); legend boxoff
l1 = legend([h8,h6,h5,h3],{'1-D Modal Theory','2-D Toy Model','1-D Toy Model','ERA5'},'Location','NorthWest'); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[-0.012 +0.02 0 0]);
l1.FontSize = 9;
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

clear('h3','h5','h6','h8')

subplot(1,2,2)

%h1 = plot(lambda_ncep2_SH_x,'Linewidth',1.5,'Color','g'); hold on
%plot(lambda_ncep2_SH_xt,'Linewidth',1.5,'Color','m','linestyle',':'); hold on;
%h2 = plot(lambda_erai_SH_x,'Linewidth',1.5,'Color','r'); hold on;
%plot(lambda_erai_SH_xt,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
h3 = plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
h5 = plot(lambda_SH_theory_1d,'linewidth',1.5,'Color','k'); hold on;
h6 = plot(lambda_SH_theory_2d,'k:','linewidth',1.5); hold on;
%h6 = plot(lambda_SH_theory_2d_refined,'k:','linewidth',1.5); hold on;
h8 = plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
%legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
%l1 = legend('NCEP2','ERAI','ERA5','1-D Toy Model','1-D Modal Theory','Location','Northwest'); legend boxoff
%l1 = legend([h7,h6,h5,h3,h4,h2,h1],{'1-D Modal Theory','1-D Full Spectrum','1-D Toy Model','ERA5','ERA5 Coarse','ERAI','NCEP2'},'Location','NorthWest'); legend boxoff
%l2 = legend([h8,h7,h6,h5,h3],{'1-D Modal Theory','2-D Toy Model Refined','2-D Toy Model','1-D Toy Model','ERA5'},'Location','NorthWest'); legend boxoff
l1 = legend([h8,h6,h5,h3],{'1-D Modal Theory','2-D Toy Model','1-D Toy Model','ERA5'},'Location','NorthWest'); legend boxoff
l1.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')