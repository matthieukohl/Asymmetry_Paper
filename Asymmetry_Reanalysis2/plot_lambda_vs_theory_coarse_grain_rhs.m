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

% era5 coarse grain

% load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p25_res.mat'); 
% load('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p25_res.mat'); 

load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p5_res.mat'); 
load('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p5_res.mat'); 

% load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p75_res.mat'); 
% load('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p75_res.mat'); 

%load('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p75_res.mat'); 

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

% non-weighted

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1));

% cosine - weighted

weight = repmat(cosd(lat(lat70:lat30)),1,1,12);

lambda_era5_NH_x_weighted = squeeze(sum(weight.*lambda_era5_NH_x_coarse_grain(lat70:lat30,p500,:),1)./sum(weight,1));
lambda_era5_NH_xy_weighted = squeeze(lambda_era5_NH_xy_coarse_grain(p500,:));

clear('level','lat','p500','lat30','lat70')

% load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p25_res.mat');
% load('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p25_res.mat');

% load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p75_res.mat');
% load('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p75_res.mat');

load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p5_res.mat');
load('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p5_res.mat');


%load('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p75_res.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

% non-weighted

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1));

% cosine - weighted

weight = repmat(cosd(lat(latm30:latm70)),1,1,12);

lambda_era5_SH_x_weighted = squeeze(sum(weight.*lambda_era5_SH_x_coarse_grain(latm30:latm70,p500,:),1)./sum(weight,1));
lambda_era5_SH_xy_weighted = squeeze(lambda_era5_SH_xy_coarse_grain(p500,:));


clear('level','lat','p500','lat30','lat70')


% load r and k values

% load('data_theory/theory_era5_NH_coarse_grained_july_5th.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th.mat','r_500_SH','k_SH_500','S_500_SH')

% load('data_theory/theory_era5_NH_coarse_grained_july_5th_res_1p75.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th_res_1p75.mat','r_500_SH','k_SH_500','S_500_SH')

% load('data_theory/theory_era5_NH_coarse_grained_july_5th_res_1p75_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_5th_res_1p75_weighted.mat','r_500_SH','k_SH_500','S_500_SH')

% load('data_theory/theory_era5_NH_coarse_grained_july_12th_res_1p25_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_july_12th_res_1p25_weighted.mat','r_500_SH','k_SH_500','S_500_SH')

% load('data_theory/theory_era5_NH_coarse_grained_res_1p5_weighted.mat','r_500_NH','k_NH_500','S_500_NH')
% load('data_theory/theory_era5_SH_coarse_grained_res_1p5_weighted.mat','r_500_SH','k_SH_500','S_500_SH')

load('data_theory/theory_era5_NH_coarse_grained_res_1p5_weighted.mat')
load('data_theory/theory_era5_SH_coarse_grained_res_1p5_weighted.mat')

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

% calculate the wavenumber from the rhs spectrum

% load lat and lon

load('lat_and_lon_era5','lon');
lon = lon(1:6:end);

[k_NH_rhs,k_SH_rhs] = deal(zeros(12,1));

R = 6371000.0; % Planetary Radius 
R_eff = cosd(45)*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);
indx = N/2+1;

% load the rhs spectra

load('spec_rhs_NH_era5_coarse_1p5_res.mat','spec_rhs_era5');
spec_NH_rhs_500 = spec_rhs_era5;
clear('spec_rhs_era5');

load('spec_rhs_SH_era5_coarse_1p5_res.mat','spec_rhs_era5');
spec_SH_rhs_500 = spec_rhs_era5;
clear('spec_rhs_era5');

for ii = 1:12

k_NH_rhs(ii) = sum(k(indx:end).*spec_NH_rhs_500(indx:end,ii)')./sum(spec_NH_rhs_500(indx:end,ii)');
k_SH_rhs(ii) = sum(k(indx:end).*spec_SH_rhs_500(indx:end,ii)')./sum(spec_SH_rhs_500(indx:end,ii)');

end

% toy-model

lambda_NH_theory = zeros(12,1);
lambda_SH_theory = zeros(12,1);

lambda_NH_theory_fixed = zeros(12,1);
lambda_SH_theory_fixed = zeros(12,1);

for ii = 1:12
    
lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_rhs(ii)*LD_NH_500(ii)); 
lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_rhs(ii)*LD_SH_500(ii));

% hold k fixed over the seasonal cycle

lambda_NH_theory_fixed(ii) = toy_model_nondimensional(r_500_NH(ii),k_NH_rhs(1)*LD_NH_500(1)); 
lambda_SH_theory_fixed(ii) = toy_model_nondimensional(r_500_SH(ii),k_SH_rhs(1)*LD_SH_500(1));

% % take k as yearly average
% 
% lambda_NH_theory_fixed(ii) = toy_model_nondimensional(r_500_NH(ii),mean(k_NH_500)*mean(LD_NH_500)); 
% lambda_SH_theory_fixed(ii) = toy_model_nondimensional(r_500_SH(ii),mean(k_SH_500)*mean(LD_SH_500));


% lambda_NH_theory(ii) = toy_model_nondimensional(r_500_NH(ii),1.8); 
% lambda_SH_theory(ii) = toy_model_nondimensional(r_500_SH(ii),1.8);


end



% modal theory

lambda_NH_modal_theory = zeros(12,1);
lambda_SH_modal_theory = zeros(12,1);

% need to look into this again; for asymmetry.mat lambda mode is not as
% smooth over the seasonal cycle as r itself. seems that there is a
% jump of lambda for certain r values when an additional 
% updraft/downdraft fits into the box; running on bigger domain helps
% but then the SH, for which r is larger, looks non-smooth

%load('mode_asymmetry_bigger_domain.mat');
%load('mode_asymmetry.mat','lambda_SH_modal_theory');

load('data_theory/theory_mode_era5_coarse_grained_res_1p5_weighted_even_bigger_domain.mat');
load('data_theory/theory_mode_era5_coarse_grained_res_1p5_weighted.mat','lambda_SH_modal_theory');


%load('mode_asymmetry.mat');


% make figure

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

h1 = plot(lambda_era5_NH_xy_weighted,'linewidth',1.5,'Color','b'); hold on;
%h11 = plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b','Linestyle','--'); hold on;

h2 = plot(lambda_NH_theory,'linewidth',1.5,'Color','k','Linestyle','-.'); hold on;
h3 = plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
l1 = legend([h3,h2,h1],{'1-D Modal Theory','1-D Toy Model','ERA5 Coarse-Grained'},'Location','NorthWest'); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[-0.01 0.02 0 0]);
l1.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

h1 = plot(lambda_era5_SH_xy_weighted,'linewidth',1.5,'Color','b'); hold on;
%h11 = plot(lambda_era5_SH_x_weighted,'linewidth',1.5,'Color','b','Linestyle','--'); hold on;
h2 = plot(lambda_SH_theory,'linewidth',1.5,'Color','k','Linestyle','-.'); hold on;
h3 = plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
l2 = legend([h3,h2,h1],{'1-D Modal Theory','1-D Toy Model','ERA5 Coarse-Grained'},'Location','NorthWest'); legend boxoff
l2pos = get(l2,'position');
set(l2,'position',l2pos+[-0.01 0.02 0 0]);
l2.FontSize = 9;
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_cycle_reanalysis_vs_theory_coarse_grain_x_mean_rhs','epsc');


figure(2);
%pos = get(hfig,'position');
%set(hfig,'position',pos.*[.5 1 2 1])
x0=10; y0=10; width=1200; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(r_500_NH,'k','Linewidth',1.5); hold on;
plot(r_500_SH,'k--','Linewidth',1.5)
legend('NH','SH','Location','SouthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Reduction factor r')
ylim([0 0.5])
yticks([0 0.1 0.2 0.3 0.4 0.5])
title('Reduction Factor r (500hPa)')
%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/r_seasonal_cycle','epsc');

subplot(1,2,2)

plot(k_NH_rhs.*LD_NH_500,'k','Linewidth',1.5); hold on;
plot(k_SH_rhs.*LD_SH_500,'k--','Linewidth',1.5)
legend('NH','SH','Location','NorthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Wavenumber k')
%yticks([0.1 0.2 0.3 0.4 0.5])
title('Wavenumber k (500hPa)')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/k_seasonal_cycle','epsc');

saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/r_and_k_seasonal_cycle_coarse_grained_rhs','epsc');

