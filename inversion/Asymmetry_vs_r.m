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

list = {'rescale0.01','rescale0.02','rescale0.05','rescale0.1','rescale0.2','rescale0.3',...                     
'rescale0.4','rescale0.6','rescale0.8','rescale1.0'};

list1 = {'Total_rhs_smooth0_iter20_day200.mat','Total_rhs_smooth0_iter20_day225.mat',...
    'Total_rhs_smooth0_iter20_day250.mat','Total_rhs_smooth0_iter30_geostrophic_w_day275.mat'};

list2 = {'RHS_rhs_smooth0_iter20_day200.mat','RHS_rhs_smooth0_iter20_day225.mat',...
    'RHS_rhs_smooth0_iter20_day250.mat','RHS_rhs_smooth0_iter20_day275.mat'};

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

strat_level_mode = 12;

lambda_mean_omega = zeros(1,length(r));
lambda_mean_omegaQG = zeros(1,length(r));
lambda_mean_omegaQG_RHS = zeros(1,length(r));

lambda_mean_omega_mode = zeros(1,length(r));
lambda_mean_omegaQG_mode = zeros(1,length(r));
lambda_mean_omegaQG_RHS_mode = zeros(1,length(r));


for tt = 1:10
tt
omega_QG = [];
omega_QG_RHS = [];
omega = [];

% Non-Linear Calculations

for i=4:4
    
s = load(strcat('output/r_simulation/',list{tt},'/',list1{i}));
t = load(strcat('output/r_simulation/',list{tt},'/',list2{i}));
omega = cat(4,omega,s.omega(:,:,:,:));
omega_QG=cat(4,omega_QG,s.omega_QG(:,:,:,:));
omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(:,:,:,:));
end

% Calculate the Asymmetry Factor Lambda

omega = omega - mean(omega(:,:,:,:),2);
lambda_omega = Lambda_weighted(omega(:,:,1:strat_level(tt),:),lat);
%lambda_omega = Skewness(omega(:,:,1:strat_level(tt),:));
%lambda_omega = Lambda_weighted(omega(:,:,6:12,:),lat);

omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
lambda_omegaQG = Lambda_weighted(omega_QG(:,:,1:strat_level(tt),:),lat);
%lambda_omegaQG = Skewness(omega_QG(:,:,1:strat_level(tt),:));
%lambda_omegaQG = Lambda_weighted(omega_QG(:,:,6:12,:),lat);

omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
lambda_omegaQG_RHS = Lambda_weighted(omega_QG_RHS(:,:,1:strat_level(tt),:),lat);
%lambda_omegaQG_RHS = Skewness(omega_QG_RHS(:,:,1:strat_level(tt),:));
%lambda_omegaQG_RHS = Lambda_weighted(omega_QG_RHS(:,:,6:12,:),lat);

lambda_mean_omega(tt) = mean(lambda_omega);
lambda_mean_omegaQG(tt) = mean(lambda_omegaQG);
lambda_mean_omegaQG_RHS(tt) = mean(lambda_omegaQG_RHS);

% Mode Calculation

% path_omega = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_',list{tt},'/','omega_day0560h00.nc');
% omega_in= ncread(path_omega,'omega_plevel', start, count);
% 
% [omega_mode] = ...
%   deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
%  
% for k = 1 : length(level_indices)
%     for t = 1 : length(event_timespan)
%         omega_mode(:, :, k, t) = omega_in(:, :, k, t)';
%     end
% end

path_mode = load(strcat('output/r_simulation_mode/',list{tt},'/Total_rhs_smooth0_iter250_day0560h00.mat'));
path_mode_RHS = load(strcat('output/r_simulation_mode/',list{tt},'/RHS_rhs_smooth0_iter20_day0560h00.mat'));
omega_mode = path_mode.omega(:,:,:,:);
omega_QG_mode = path_mode.omega_QG(:,:,:,:);
omega_QG_RHS_mode = path_mode_RHS.omega_QG(:,:,:,191:194);

omega_mode = omega_mode - mean(omega_mode(:,:,:,:),2);
lambda_omega_mode = Lambda_weighted(omega_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omega_mode = Skewness(omega_mode(:,:,1:strat_level_mode,:));
%lambda_omega_mode = Lambda_weighted(omega_mode(:,:,6:12,:),lat);

omega_QG_mode = omega_QG_mode-mean(omega_QG_mode(:,:,:,:),2);
lambda_omegaQG_mode = Lambda_weighted(omega_QG_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omegaQG_mode = Skewness(omega_QG_mode(:,:,1:strat_level_mode,:));
%lambda_omegaQG_mode = Lambda_weighted(omega_QG_mode(:,:,6:12,:),lat);

omega_QG_RHS_mode = omega_QG_RHS_mode-mean(omega_QG_RHS_mode(:,:,:,:),2);
lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,1:strat_level_mode,:),lat);
%lambda_omegaQG_RHS_mode = Skewness(omega_QG_RHS_mode(:,:,1:strat_level_mode,:));
%lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,6:12,:),lat);

lambda_mean_omega_mode(tt) = mean(lambda_omega_mode);
lambda_mean_omegaQG_mode(tt) = mean(lambda_omegaQG_mode);
lambda_mean_omegaQG_RHS_mode(tt) = mean(lambda_omegaQG_RHS_mode);

end

lambda_theory_mode = [0.5002, 0.5094, 0.5271, 0.5452, 0.5709, 0.6003, 0.6372, 0.6810,...
0.7448, 0.7547, 0.7733, 0.7882, 0.8063, 0.8169, 0.8296, 0.8520, 0.8755, 0.9007,0.9317, 0.9530];

%lambda_theory_mode = [0.5002,0.5155,0.5278,0.5526,0.5793,0.6132,0.6543,0.7012...
 %0.7763, 0.7917, 0.8067, 0.8244, 0.8413, 0.8508, 0.8598, 0.8772, 0.8957, 0.9160,0.9405,0.9578];

r_theory = [1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.18,0.16,0.14,0.12,0.11,0.10, ...
    0.08,0.06,0.04,0.02,0.01];

% lambda_theory_turbulent = [0.7084, 0.6734, 0.6529, 0.6243, 0.5870, 0.5632, ...
%  0.5482, 0.5281, 0.5118, 0.5000];

%lambda_theory_turbulent = [0.7612,0.7431,0.7082,0.6727,0.6289,0.5997,0.5775,0.5443,0.5196,0.5000];

lambda_theory_turbulent_randn = [0.7585,0.7585,0.7578,0.6941,0.6858,0.6528,0.5987,0.5438,0.5177,0.4890];

lambda_theory_turbulent_randn_30 = [0.7535, 0.7535,0.7518,0.7184,0.6639,0.6259,...  
0.6038, 0.5383, 0.5217,0.5003];

lambda_theory_turbulent_randn_50 = [0.766,0.7666,0.7420,0.7063,0.6668,0.6340,0.6183,0.5595,0.5313,0.5164];

figure(1)
semilogx(r,lambda_mean_omega); hold on;
semilogx(r,lambda_mean_omegaQG); hold on;
semilogx(r,lambda_mean_omegaQG_RHS); hold on
%semilogx(r,lambda_mean_omega_mode); hold on;
%semilogx(r,lambda_mean_omegaQG_mode); hold on;
%semilogx(r,lambda_mean_omegaQG_RHS_mode); hold on;
%semilogx(r_theory,lambda_theory_mode); hold on;
semilogx(r,lambda_theory_turbulent_randn_50);
xlabel('Optical Thickness')
ylabel('Skewness')
title('Asymmetry with w lower boundary')
legend('Omega','QG','QG_{RHS}','Omega_{mode}','QG_{mode}','QG_{RHS_Mode}','Theory_{Mode}','Theory_{Turbulent}');

savefig('Figure_Asymmetry_r_theory5.fig');

 