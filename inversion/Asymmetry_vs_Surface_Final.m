clear;
close all;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output_abs/abs1.0_T85/temp_day1450h00.nc';
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);

lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;
plevels = double(ncread(path_temp, 'level'));
level_indices = [5:30];
level = plevels(level_indices);


event_latspan = [83:111]; 
event_lonspan = [1:256];
event_timespan = [1:200];

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

list1 = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

abs = [0.4,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,4.0,6.0];

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

strat_level_mode =  [8; 9; 10; 10 ; 11; 11; 12; 12; 13; 13; 14; 15; 15; 16; 18];

%strat_level_mode = [12,13,14,14,15,15,16,16,17,17,18,19,19,20,22];

strat_level_turbulent = [8,9,10,10,11,11,12,12,13,13,14,14,15,16,18];

lambda_mean_omega = zeros(1,15);
lambda_mean_omegaQG = zeros(1,15);
lambda_mean_omegaQG_RHS = zeros(1,15);
lambda_mean_omega_mode = zeros(1,15);
lambda_mean_omegaQG_mode = zeros(1,15);
lambda_mean_omegaQG_RHS_mode = zeros(1,15);

for tt = 1:15
tt
omega_QG = [];
omega_QG_RHS = [];
omega = [];
r = [];

% % Load Non-Linear Simulations True Boundary
% 
% for i=1:6
%     
% s = load(strcat('output/abs_simulation/',list{tt},'/final/Total_rhs_smooth0_rms4_geostrophic_w_trial_day',num2str(5 *(28+i)),'.mat'));
% t = load(strcat('output/abs_simulation/',list1{tt},'/final/RHS_rhs_smooth0_rms4_geostrophic_w_trial_day',num2str(5 *(28+i)),'.mat'));
% omega_QG=cat(4,omega_QG,s.omega_QG(11:26,:,:,:));
% omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(11:26,:,:,:));
% omega = cat(4,omega,s.omega(11:26,:,:,:));
% r = cat(4,r,t.r);
% 
% end

% Load Non-Linear Simulations Zero Boundary

for i=1:6
    
s = load(strcat('output/abs_simulation/',list{tt},'/final/Total_rhs_smooth0_rms4_geostrophic_trial_day',num2str(5 *(28+i)),'.mat'));
t = load(strcat('output/abs_simulation/',list1{tt},'/final/RHS_rhs_smooth0_rms4_geostrophic_trial_day',num2str(5 *(28+i)),'.mat'));
omega_QG=cat(4,omega_QG,s.omega_QG(11:26,:,:,:));
omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(11:26,:,:,:));
omega = cat(4,omega,s.omega(11:26,:,:,:));
r = cat(4,r,t.r);

end

% Load Mode Simulations: True BC

%path_mode = load(strcat('output/abs_mode_simulation/',list{tt},'/final/Total_rhs_smooth0_rms15_geostrophic_w_trial_day0560h00.mat'));
%path_mode_RHS = load(strcat('output/abs_mode_simulation/',list{tt},'/final/RHS_rhs_smooth0_rms15_geostrophic_w_trial_day0560h00.mat'));
%path_mode_LHS = load(strcat('output/abs_mode_simulation/',list{tt},'/final/LHS_rhs_smooth0_rms4_geostrophic_day0560h00.mat'))
%path_mode = load(strcat('output/abs_mode_simulation/',list{tt},'/Inversion_smooth0_it30_geostrophic_day.mat'));
%path_mode_RHS = load(strcat('output/abs_mode_simulation/',list{tt},'/Inversion_RHS_smooth0_it30_geostrophic_w1_day.mat'));

% Load Mode Simulations: Zero BC

path_mode = load(strcat('output/abs_mode_simulation/',list{tt},'/final/Total_rhs_smooth0_rms15_geostrophic_trial_day0560h00.mat'));
path_mode_RHS = load(strcat('output/abs_mode_simulation/',list{tt},'/final/RHS_rhs_smooth0_rms15_geostrophic_trial_day0560h00.mat'));



omega_mode = path_mode.omega(11:26,:,:,:);
omegaQG_mode = path_mode.omega_QG(11:26,:,:,:);
omegaQG_RHS_mode = path_mode_RHS.omega_QG(11:26,:,:,:);

% % Calculate the Asymmetry Factor Lambda
% 
omega = omega - mean(omega(:,:,:,:),2);
lambda_omega = Lambda(omega(:,:,1:strat_level_turbulent(tt),:));
% %lambda_omega = Skewness(omega(:,:,1:strat_level_turbulent(tt),:));
%lambda_omega = Lambda(omega(:,:,6,:));
% 
omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
lambda_omegaQG = Lambda(omega_QG(:,:,1:strat_level_turbulent(tt),:));
% %lambda_omegaQG = Skewness(omega_QG(:,:,1:strat_level_turbulent(tt),:));
%lambda_omegaQG = Lambda(omega_QG(:,:,6,:));
% 
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS(:,:,1:strat_level_turbulent(tt),:));
% %lambda_omegaQG_RHS = Skewness(omega_QG_RHS(:,:,1:strat_level_turbulent(tt),:));
%lambda_omegaQG_RHS = Lambda(omega_QG_RHS(:,:,6,:));

omega_mode = omega_mode-mean(omega_mode(:,:,:,:),2);
lambda_omega_mode = Lambda(omega_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omega_mode = Skewness(omega_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omega_mode = Lambda(omega_mode(:,:,6,:));

omegaQG_mode = omegaQG_mode - mean(omegaQG_mode(:,:,:,:),2);
lambda_omegaQG_mode = Lambda(omegaQG_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omegaQG_mode = Skewness(omegaQG_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omegaQG_mode = Lambda(omegaQG_mode(:,:,6,:));

omegaQG_RHS_mode = omegaQG_RHS_mode - mean(omegaQG_RHS_mode(:,:,:,:),2);
lambda_omegaQG_RHS_mode = Lambda(omegaQG_RHS_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omegaQG_RHS_mode = Skewness(omegaQG_RHS_mode(:,:,1:strat_level_mode(tt),:));
%lambda_omegaQG_RHS_mode = Lambda(omegaQG_RHS_mode(:,:,6,:));


lambda_mean_omega(tt) = mean(lambda_omega);
lambda_mean_omegaQG(tt) = mean(lambda_omegaQG);
lambda_mean_omegaQG_RHS(tt) = mean(lambda_omegaQG_RHS);

lambda_mean_omega_mode(tt) = mean(lambda_omega_mode);
lambda_mean_omegaQG_mode(tt) = mean(lambda_omegaQG_mode);
lambda_mean_omegaQG_RHS_mode(tt) = mean(lambda_omegaQG_RHS_mode);

end

% lambda_theory = [0.5149, 0.5252, 0.5312, 0.5367, 0.5432, 0.5509, 0.5646, ...
% 0.5785, 0.5935, 0.6090, 0.6197, 0.6453, 0.6668, 0.6809, 0.6868];

% lambda_theory_rand = [0.5043,0.5486,0.5420,0.5718,0.4612,0.6157,0.5167,0.5223 ...    
% 0.6337, 0.6420, 0.6818, 0.5922, 0.6497, 0.6451, 0.8244];

lambda_theory_rand = [0.4710,0.5560,0.5690,0.5626,0.4330,0.6503,0.5888,...
0.6297,0.6351,0.6267,0.6415,0.7451,0.6629,0.7939, 0.8061];



% figure(1)
% plot(surf_temp,lambda_mean_omega,'Color','blue'); hold on;
% plot(surf_temp,lambda_mean_omegaQG,'Color','blue','LineStyle','--'); hold on;
% plot(surf_temp,lambda_mean_omegaQG_RHS,'Color','blue','LineStyle','-.'); hold on;
% % %plot(surf_temp,lambda_theory_rand);
% xlabel('Surface Temperature')
% ylabel('Lambda')
% title('Asymmetry Geostrophic')
% legend('Omega','QG','QG_{RHS}','Theory');
% 
% 
% figure(2)
% plot(surf_temp,lambda_mean_omega_mode,'Color','red'); hold on;
% plot(surf_temp,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--'); hold on;
% plot(surf_temp,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','-.'); hold on;
% xlabel('Surface Temperature')
% ylabel('Lambda')
% title('Asymmetry Geostrophic')
% legend('Omega_{mode}','QG_{mode}','QG_{RHS_Mode}')
%savefig('Figure_Asymmetry_vs_Surface_Temperature_Mode1.fig');

% save('Asymmetry_vs_abs_zero_data.mat','lambda_mean_omega','lambda_mean_omegaQG',...
%     'lambda_mean_omegaQG_RHS','lambda_mean_omega_mode','lambda_mean_omegaQG_mode',...
%     'lambda_mean_omegaQG_RHS_mode')
% 

file = load('Asymmetry_vs_abs_data.mat');
lambda_mean_omega = file.lambda_mean_omega;
lambda_mean_omegaQG = file.lambda_mean_omegaQG;
lambda_mean_omegaQG_RHS = file.lambda_mean_omegaQG_RHS;
lambda_mean_omega_mode = file.lambda_mean_omega_mode;
lambda_mean_omegaQG_mode = file.lambda_mean_omegaQG_mode;
lambda_mean_omegaQG_RHS_mode = file.lambda_mean_omegaQG_RHS_mode;

figure(3)
plot(surf_temp,lambda_mean_omega_mode,'Color','red','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
%plot(surf_temp,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omega,'Color','blue','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omegaQG,'Color','blue','LineStyle','--','linewidth',1.2); hold on;
%plot(surf_temp,lambda_mean_omegaQG_RHS,'Color','blue','Linestyle','--','linewidth',1.2); hold on;
xlabel('Global-mean surface air temperature (K)')
ylabel('Asymmetry parameter \lambda')
%legend('Omega','QG','QG_{RHS}','Omega_{mode}','QG_{mode}','QG_{RHS_{Mode}}','Location','northwest')
%legend('GCM_{Mode}','QG_{Mode}(RHS)','GCM_{Turbul.}','QG_{Turbul.}(RHS)','Location','northwest')
legend('GCM_{Mode}','QG_{Mode}','GCM_{Turbul.}','QG_{Turbul.}','Location','northwest')
%legend('Omega','QG_{RHS}','Omega_{mode}','QG_{RHS_{Mode}}','Location','northwest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
ylim([0.4,1])
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Asymmetry_vs_abs_full_zero','epsc')
close all;

figure(4)
plot(surf_temp,lambda_mean_omega_mode,'Color','red','linewidth',1.2); hold on;
%plot(surf_temp,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omega,'Color','blue','linewidth',1.2); hold on;
%plot(surf_temp,lambda_mean_omegaQG,'Color','blue','LineStyle','--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mean_omegaQG_RHS,'Color','blue','Linestyle','--','linewidth',1.2); hold on;
xlabel('Global-mean surface air temperature (K)')
ylabel('Asymmetry parameter \lambda')
%legend('Omega','QG','QG_{RHS}','Omega_{mode}','QG_{mode}','QG_{RHS_{Mode}}','Location','northwest')
legend('GCM_{Mode}','QG_{Mode}(RHS)','GCM_{Turbul.}','QG_{Turbul.}(RHS)','Location','northwest')
%legend('GCM_{Mode}','QG_{Mode}','GCM_{Turbul.}','QG_{Turbul.}','Location','northwest')
%legend('Omega','QG_{RHS}','Omega_{mode}','QG_{RHS_{Mode}}','Location','northwest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
ylim([0.4,1])
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Asymmetry_vs_abs_RHS_zero','epsc')
close all;

% save('Asymmetry_vs_abs_data.mat','lambda_mean_omega','lambda_mean_omegaQG',...
% 'lambda_mean_omegaQG_RHS','lambda_mean_omega_mode','lambda_mean_omegaQG_mode',...
% 'lambda_mean_omegaQG_RHS_mode')

%print('Plots_Generals/test','-dpdf')


