clear;
close all;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale1.0/temp_day0025h00.nc';
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
level_indices = [1:30];
level = plevels(level_indices);


event_latspan = [83:111]; 
event_lonspan = [1:256];
event_timespan = [1:200];

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;

list = {'rescale0.01','rescale0.02','rescale0.05','rescale0.1','rescale0.2','rescale0.3',...                     
'rescale0.4','rescale0.6','rescale0.8','rescale1.0'};

list_t = {'rescale_true0.01','rescale0.02','rescale0.05','rescale0.1','rescale0.2','rescale0.3',...                     
'rescale0.4','rescale_true0.6','rescale0.8','rescale_true1.0'};

%list1 = {'Total_rhs_smooth0_rms4_geostrophic_w_day225.mat','Total_rhs_smooth0_rms4_geostrophic_w_day250.mat',...
    %'Total_rhs_smooth0_rms4_geostrophic_w_day275.mat','Total_rhs_smooth0_rms4_geostrophic_w_day300.mat'};

%list2 = {'RHS_rhs_smooth0_rms4_geostrophic_w_day225.mat','RHS_rhs_smooth0_rms4_geostrophic_w_day250.mat',...
    %'RHS_rhs_smooth0_rms4_geostrophic_w_day275.mat','RHS_rhs_smooth0_rms4_geostrophic_w_day300.mat'};

list1 = {'Total_rhs_smooth0_rms4_geostrophic_day225.mat','Total_rhs_smooth0_rms4_geostrophic_day250.mat',...
    'Total_rhs_smooth0_rms4_geostrophic_day275.mat','Total_rhs_smooth0_rms4_geostrophic_day300.mat'};

list2 = {'RHS_rhs_smooth0_rms4_geostrophic_day225.mat','RHS_rhs_smooth0_rms4_geostrophic_day250.mat',...
    'RHS_rhs_smooth0_rms4_geostrophic_day275.mat','RHS_rhs_smooth0_rms4_geostrophic_day300.mat'};




r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

strat_level_RHS = 4+[19;14;13;14;14;13;12;18;12;18];

strat_level = [18,18,18,18,17,17,16,16,16,16];

strat_level_mode = 12+4;


lambda_mean_omega = zeros(1,length(r));
lambda_mean_omegaQG = zeros(1,length(r));
lambda_mean_omegaQG_RHS = zeros(1,length(r));

lambda_mean_omega_mode = zeros(1,length(r));
lambda_mean_omegaQG_mode = zeros(1,length(r));
lambda_mean_omegaQG_RHS_mode = zeros(1,length(r));
lambda_mean_omegaQG_LHS_mode = zeros(1,length(r));


for tt = 1:10
tt
omega_QG = [];
omega_QG_RHS = [];
omega = [];

%Non-Linear Calculations

for i=1:4
    
s = load(strcat('output/r_simulation/',list_t{tt},'/final/',list1{i}));
t = load(strcat('output/r_simulation/',list{tt},'/final/',list2{i}));
omega = cat(4,omega,s.omega(:,:,:,:));
omega_QG=cat(4,omega_QG,s.omega_QG(:,:,:,:));
omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(:,:,:,:));
end

% omega = omega(:,:,5:30,:);
% omega_QG = omega_QG(:,:,5:30,:);
% omega_QG_RHS = omega_QG_RHS(:,:,5:30,:);

% Calculate the Asymmetry Factor Lambda

omega = omega - mean(omega(:,:,:,:),2);
%lambda_omega = Lambda_weighted(omega(:,:,1:strat_level(tt),:),lat);
%lambda_omega = Skewness(omega(:,:,1:strat_level(tt),:));
lambda_omega = Lambda_weighted(omega(:,:,11,:),lat);

%omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
%lambda_omegaQG = Lambda_weighted(omega_QG(:,:,1:strat_level(tt),:),lat);
%lambda_omegaQG = Skewness(omega_QG(:,:,1:strat_level(tt),:));
lambda_omegaQG = Lambda_weighted(omega_QG(:,:,11,:),lat);
%lambda_omegaQG = Lambda(omega_QG(:,:,1:strat_level(tt),:),lat);

%omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
%lambda_omegaQG_RHS = Lambda_weighted(omega_QG_RHS(:,:,1:strat_level_RHS(tt),:),lat);
%lambda_omegaQG_RHS = Skewness(omega_QG_RHS(:,:,1:strat_level(tt),:));
lambda_omegaQG_RHS = Lambda_weighted(omega_QG_RHS(:,:,11,:),lat);
%lambda_omegaQG_RHS = Lambda(omega_QG_RHS(:,:,1:strat_level_RHS(tt),:),lat);

lambda_mean_omega(tt) = mean(lambda_omega);
lambda_mean_omegaQG(tt) = mean(lambda_omegaQG);
lambda_mean_omegaQG_RHS(tt) = mean(lambda_omegaQG_RHS);

% Mode Calculation

%True-Boundary
%path_mode = load(strcat('output/r_simulation_mode/',list{tt},'/final/Total_rhs_smooth0_rms6_geostrophic_w_day0560h00.mat'));
%path_mode_RHS = load(strcat('output/r_simulation_mode/',list{tt},'/final/RHS_rhs_smooth0_rms6_geostrophic_w_day0560h00.mat'));
%Zero-Boundary
path_mode = load(strcat('output/r_simulation_mode/',list{tt},'/final/Total_rhs_smooth0_rms4_geostrophic_day0560h00.mat'));
path_mode_RHS = load(strcat('output/r_simulation_mode/',list{tt},'/final/RHS_rhs_smooth0_rms4_geostrophic_day0560h00.mat'));

%path_mode_LHS = load(strcat('output/r_simulation_mode/',list{tt},'/final/LHS_rhs_smooth0_rms6_geostrophic_w_day0560h00.mat'));

omega_mode = path_mode.omega(:,:,:,:);
omega_QG_mode = path_mode.omega_QG(:,:,:,:);
omega_QG_RHS_mode = path_mode_RHS.omega_QG(:,:,:,:);
%omega_QG_LHS_mode = path_mode_LHS.omega_QG(:,:,:,:);

% omega_mode = omega_mode(:,:,5:30,:);
% omega_QG_mode = omega_QG_mode(:,:,5:30,:);
% omega_QG_RHS_mode = omega_QG_RHS_mode(:,:,5:30,:);

%omega_mode = omega_mode - mean(omega_mode(:,:,:,:),2);
%lambda_omega_mode = Lambda_weighted(omega_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omega_mode = Skewness(omega_mode(:,:,1:strat_level_mode,:));
lambda_omega_mode = Lambda_weighted(omega_mode(:,:,11,:),lat);
%lambda_omega_mode = Lambda(omega_mode(:,:,1:strat_level_mode,:),lat);

%omega_QG_mode = omega_QG_mode-mean(omega_QG_mode(:,:,:,:),2);
%lambda_omegaQG_mode = Lambda_weighted(omega_QG_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omegaQG_mode = Skewness(omega_QG_mode(:,:,1:strat_level_mode,:));
lambda_omegaQG_mode = Lambda_weighted(omega_QG_mode(:,:,11,:),lat);
%lambda_omegaQG_mode = Lambda(omega_QG_mode(:,:,1:strat_level_mode,:),lat);

%omega_QG_RHS_mode = omega_QG_RHS_mode-mean(omega_QG_RHS_mode(:,:,:,:),2);
%lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omegaQG_RHS_mode = Skewness(omega_QG_RHS_mode(:,:,1:strat_level_mode,:));
lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,11,:),lat);
%lambda_omegaQG_RHS_mode = Lambda(omega_QG_RHS_mode(:,:,1:strat_level_mode,:),lat);

%omega_QG_LHS_mode = omega_QG_LHS_mode-mean(omega_QG_LHS_mode(:,:,:,:),2);
%lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omegaQG_RHS_mode = Skewness(omega_QG_RHS_mode(:,:,1:strat_level_mode,:));
%lambda_omegaQG_LHS_mode = Lambda_weighted(omega_QG_LHS_mode(:,:,11,:),lat);

lambda_mean_omega_mode(tt) = mean(lambda_omega_mode);
lambda_mean_omegaQG_mode(tt) = mean(lambda_omegaQG_mode);
lambda_mean_omegaQG_RHS_mode(tt) = mean(lambda_omegaQG_RHS_mode);
%lambda_mean_omegaQG_LHS_mode(tt) = mean(lambda_omegaQG_LHS_mode);

end

% lambda_theory_mode = [0.5002, 0.5094, 0.5271, 0.5452, 0.5709, 0.6003, 0.6372, 0.6810,...
% 0.7448, 0.7547, 0.7733, 0.7882, 0.8063, 0.8169, 0.8296, 0.8520, 0.8755, 0.9007,0.9317, 0.9530];

%lambda_theory_mode = [0.5002,0.5155,0.5278,0.5526,0.5793,0.6132,0.6543,0.7012...
 %0.7763, 0.7917, 0.8067, 0.8244, 0.8413, 0.8508, 0.8598, 0.8772, 0.8957, 0.9160,0.9405,0.9578];
 
 lambda_theory_mode_8_N400 = [0.4994, 0.5168, 0.5371, 0.5573, 0.5895, 0.6326, 0.6630, 0.6957, 0.7629,...
  0.7765,0.7952,0.8127,0.8303,0.8397,0.8218,0.8450,0.8688,0.8960,0.9302,0.9302];

r_theory = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.18,0.16,0.14,0.12,0.11,0.10, ...
    0.08,0.06,0.04,0.02,0.01];

% lambda_theory_turbulent = [0.7084, 0.6734, 0.6529, 0.6243, 0.5870, 0.5632, ...
%  0.5482, 0.5281, 0.5118, 0.5000];

%lambda_theory_turbulent = [0.7612,0.7431,0.7082,0.6727,0.6289,0.5997,0.5775,0.5443,0.5196,0.5000];

lambda_theory_turbulent_randn = [0.7585,0.7585,0.7578,0.6941,0.6858,0.6528,0.5987,0.5438,0.5177,0.4890];

lambda_theory_turbulent_randn_30 = [0.7535, 0.7535,0.7518,0.7184,0.6639,0.6259,...  
0.6038, 0.5383, 0.5217,0.5003];

lambda_theory_turbulent_randn_50 = [0.766,0.7666,0.7420,0.7063,0.6668,0.6340,0.6183,0.5595,0.5313,0.5164];

lambda_theory_correct_8_N100 = [0.4965,0.4968,0.5366,0.5369,0.5467,0.5729,...
0.6383,0.7351,0.7955,0.8246,0.8248,0.8408,0.8585,0.8658,0.8713,0.8906,...
0.9108,0.9261,0.9509,0.9732];

lambda_theory_correct_N200 = [0.5018,0.5171,0.5352,0.5504,0.5825,0.6114,...
0.6502,0.7085,0.7790,0.7883,0.8112,0.8265,0.8433,0.8512,0.8597,0.8762,...
0.8956,0.9167,0.9415,0.9590];

% file = load('Asymmetry_vs_r_data.mat');
% lambda_mean_omega = file.lambda_mean_omega;
% lambda_mean_omegaQG = file.lambda_mean_omegaQG;
% lambda_mean_omegaQG_RHS = file.lambda_mean_omegaQG_RHS;
% lambda_mean_omega_mode = file.lambda_mean_omega_mode;
% lambda_mean_omegaQG_mode = file.lambda_mean_omegaQG_mode;
% lambda_mean_omegaQG_RHS_mode = file.lambda_mean_omegaQG_RHS_mode;

figure(1)
semilogx(r,lambda_mean_omega_mode,'Color','red','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
%semilogx(r,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omega,'Color','blue','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omegaQG,'Color','blue','LineStyle','--','linewidth',1.2); hold on;
%semilogx(r,lambda_mean_omegaQG_RHS,'Color','blue','LineStyle','--','linewidth',1.2); hold on
%semilogx(r_theory,lambda_theory_correct_N200); hold on;
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
%legend('Omega','QG','QG_{RHS}','Omega_{mode}','QG_{mode}','QG_{RHS_{Mode}}');
%legend('GCM','QG','QG_{r=1}');
%legend('GCM_{mode}','QG_{mode}(RHS)','GCM_{Turbul.}','QG_{Turbul.}(RHS)');
legend('GCM_{mode}','QG_{mode}','GCM_{Turbul.}','QG_{Turbul.}');
%legend('Omega','QG_{RHS}','Omega_{mode}','QG_{RHS_{Mode}}');
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5)
ylim([0.4,1])
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/lambda_vs_r_full_zero','epsc')

figure(2)
semilogx(r,lambda_mean_omega_mode,'Color','red','linewidth',1.2); hold on;
%semilogx(r,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','--','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omega,'Color','blue','linewidth',1.2); hold on;
%semilogx(r,lambda_mean_omegaQG,'Color','blue','LineStyle','--','linewidth',1.2); hold on;
semilogx(r,lambda_mean_omegaQG_RHS,'Color','blue','LineStyle','--','linewidth',1.2); hold on
%semilogx(r_theory,lambda_theory_correct_N200); hold on;
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
%legend('Omega','QG','QG_{RHS}','Omega_{mode}','QG_{mode}','QG_{RHS_{Mode}}');
%legend('GCM','QG','QG_{r=1}');
legend('GCM_{mode}','QG_{mode}(RHS)','GCM_{Turbul.}','QG_{Turbul.}(RHS)','Location','NorthEast');
%legend('GCM_{mode}','QG_{mode}','GCM_{Turbul.}','QG_{Turbul.}');
%legend('Omega','QG_{RHS}','Omega_{mode}','QG_{RHS_{Mode}}');
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
set(gca,'linewidth',1.5);
%ylim([0.5,1]);
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/lambda_vs_r_RHS_zero','epsc')


% save('Asymmetry_vs_r_data.mat','lambda_mean_omega','lambda_mean_omegaQG',...
%    'lambda_mean_omegaQG_RHS','lambda_mean_omega_mode','lambda_mean_omegaQG_mode',...
%     'lambda_mean_omegaQG_RHS_mode')



% figure(1); saveas(gcf,'Plots_Generals/Asymmetry_vs_r/Asymmetry_vs_r_w_0.fig')
% 
% 
% figure(2)
% semilogx(r,lambda_mean_omega,'Color','blue'); hold on;
% semilogx(r,lambda_mean_omegaQG,'Color','blue','LineStyle','--'); hold on;
% semilogx(r,lambda_mean_omegaQG_RHS,'Color','blue','LineStyle','-.'); hold on
% xlabel('r')
% ylabel('Lambda')
% title('Asymmetry r-Simulation Turbulent rms6 all days')
% legend('Omega','QG','QG_{RHS}');

%figure(3)
%semilogx(r,lambda_mean_omega_mode,'Color','red'); hold on;
%semilogx(r,lambda_mean_omegaQG_mode,'Color','red','LineStyle','--'); hold on;
%semilogx(r,lambda_mean_omegaQG_RHS_mode,'Color','red','LineStyle','-.'); hold on;
%semilogx(r,lambda_mean_omegaQG_LHS_mode,'Color','red','LineStyle',':'); hold on;
%xlabel('r')
%ylabel('Lambda')
%title('Asymmetry r-Simulation Mode rms6 all days')
%legend('Omega_{mode}','QG_{mode}','QG_{RHS_Mode}','QG_{LHS_Mode}');


% Data