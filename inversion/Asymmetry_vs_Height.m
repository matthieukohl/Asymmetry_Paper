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
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808 ...
310.7297, 316.1383];

strat_level_mode =  [8; 9; 10; 10 ; 11; 11; 12; 12; 13; 13; 14; 15; 15; 16; 18];

strat_level_turbulent = [8,9,10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15,16,18];

lambda_mean_omega = zeros(1,15);
lambda_mean_omegaQG = zeros(1,15);
lambda_mean_omegaQG_RHS = zeros(1,15);
lambda_mean_omega_mode = zeros(1,15);
lambda_mean_omegaQG_mode = zeros(1,15);
lambda_mean_omegaQG_RHS_mode = zeros(1,15);

for tt = 15:15

omega_QG = [];
omega_QG_RHS = [];
omega = [];
r = [];

% Load Non-Linear Simulations

for i=1:6
    
s = load(strcat('output/abs_simulation/',list{tt},'/Inversion_Total_smooth0_day',num2str(5 *(28+i)),'.mat'));
t = load(strcat('output/abs_simulation/',list1{tt},'/Inversion_RHS_smooth0_day',num2str(5 *(28+i)),'.mat'));
omega_QG=cat(4,omega_QG,s.omega_QG(11:26,:,:,:));
omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(11:26,:,:,:));
omega = cat(4,omega,s.omega(11:26,:,:,:));
r = cat(4,r,t.r);

end

figure(1)
contour(r(:,:,end,end)); colorbar


% Load Mode Simulations

path_mode = load(strcat('output/abs_mode_simulation/',list{tt},'/Inversion_smooth0_it30_geostrophic_day.mat'));
path_mode_RHS = load(strcat('output/abs_mode_simulation/',list{tt},'/Inversion_RHS_smooth5_day.mat'));
omega_mode = path_mode.omega(11:26,:,:,:);
omegaQG_mode = path_mode.omega_QG(11:26,:,:,:);
omegaQG_RHS_mode = path_mode_RHS.omega_QG(11:26,:,:,6:9);
r_mode = path_mode.r(11:26,:,:,:);

% Calculate the Asymmetry Factor Lambda

omega = omega - mean(omega(:,:,:,:),2);
%lambda_omega = Lambda(omega(:,:,1:strat_level_turbulent(tt),:));
lambda_omega = Lambda_z(omega(:,:,:,:));

omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
%lambda_omegaQG = Lambda(omega_QG(:,:,1:strat_level_turbulent(tt),:));
lambda_omegaQG = Lambda_z(omega_QG(:,:,:,:));

omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
%lambda_omegaQG_RHS = Lambda(omega_QG_RHS(:,:,1:strat_level_turbulent(tt),:));
lambda_omegaQG_RHS = Lambda_z(omega_QG_RHS(:,:,:,:));

omega_mode = omega_mode-mean(omega_mode(:,:,:,:),2);
%lambda_omega_mode = Lambda(omega_mode(:,:,1:strat_level_mode(tt),:));
lambda_omega_mode = Lambda_z(omega_mode(:,:,:,:));

omegaQG_mode = omegaQG_mode - mean(omegaQG_mode(:,:,:,:),2);
%lambda_omegaQG_mode = Lambda(omegaQG_mode(:,:,1:strat_level_mode(tt),:));
lambda_omegaQG_mode = Lambda_z(omegaQG_mode(:,:,:,:));

omegaQG_RHS_mode = omegaQG_RHS_mode - mean(omegaQG_RHS_mode(:,:,:,:),2);
%lambda_omegaQG_RHS_mode = Lambda(omegaQG_RHS_mode(:,:,1:strat_level_mode(tt),:));
lambda_omegaQG_RHS_mode = Lambda_z(omegaQG_RHS_mode(:,:,:,:));

r = squeeze(mean(mean(mean(r,1),2),4));
r_mode = squeeze(mean(mean(mean(r_mode,1),2),4));

figure(2)
plot(lambda_omega_mode,level); hold on;
plot(lambda_omegaQG_mode,level); hold on;
plot(lambda_omegaQG_RHS_mode,level); hold on;
plot(r_mode,level)
legend('Omega','QG','RHS')
title('Asymmetry and r factor for abs=2.0 Mode')
set(gca, 'YDir','reverse')

figure(3)
plot(lambda_omega,level); hold on;
plot(lambda_omegaQG,level); hold on;
plot(lambda_omegaQG_RHS,level); hold on;
plot(r,level)
legend('Omega','QG','RHS','r')
title('Asymmetry and r factor for abs=2.0 Turbulent')
set(gca, 'YDir','reverse')


end


