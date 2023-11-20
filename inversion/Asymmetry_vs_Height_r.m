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
    'Total_rhs_smooth0_iter20_day250.mat','Total_rhs_smooth0_iter20_day275.mat'};

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

for tt = 4:4
tt
omega_QG = [];
omega_QG_RHS = [];
omega = [];

% Non-Linear Calculations

for i=1:4
    
s = load(strcat('output/r_simulation/',list{tt},'/',list1{i}));
t = load(strcat('output/r_simulation/',list{tt},'/',list2{i}));
omega = cat(4,omega,s.omega(:,:,:,:));
omega_QG=cat(4,omega_QG,s.omega_QG(:,:,:,:));
omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG(:,:,:,:));
end

% Calculate the Asymmetry Factor Lambda

omega = omega - mean(omega(:,:,:,:),2);
%lambda_omega = Lambda_weighted(omega(:,:,1:strat_level(tt),:),lat);
lambda_omega = Lambda_z(omega(:,:,:,:));

omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
%lambda_omegaQG = Lambda_weighted(omega_QG(:,:,1:strat_level(tt),:),lat);
lambda_omegaQG = Lambda_z(omega_QG(:,:,:,:));

omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
%lambda_omegaQG_RHS = Lambda_weighted(omega_QG_RHS(:,:,1:strat_level(tt),:),lat);
lambda_omegaQG_RHS = Lambda_z(omega_QG_RHS(:,:,:,:));

% Mode Calculation

path_mode = load(strcat('output/r_simulation_mode/',list{tt},'/Total_rhs_smooth0_iter250_day0560h00.mat'));
path_mode_RHS = load(strcat('output/r_simulation_mode/',list{tt},'/RHS_rhs_smooth0_iter20_day0560h00.mat'));
omega_mode = path_mode.omega(:,:,:,:);
omega_QG_mode = path_mode.omega_QG(:,:,:,:);
omega_QG_RHS_mode = path_mode_RHS.omega_QG(:,:,:,191:194);

omega_mode = omega_mode - mean(omega_mode(:,:,:,:),2);
%lambda_omega_mode = Lambda_weighted(omega_mode(:,:,1:strat_level_mode,:),lat);
lambda_omega_mode = Lambda_z(omega_mode(:,:,:,:));

omega_QG_mode = omega_QG_mode-mean(omega_QG_mode(:,:,:,:),2);
%lambda_omegaQG_mode = Lambda_weighted(omega_QG_mode(:,:,1:strat_level_mode,:),lat);
lambda_omegaQG_mode = Lambda_z(omega_QG_mode(:,:,:,:));

omega_QG_RHS_mode = omega_QG_RHS_mode-mean(omega_QG_RHS_mode(:,:,:,:),2);
%lambda_omegaQG_RHS_mode = Lambda_weighted(omega_QG_RHS_mode(:,:,1:strat_level_mode,:),lat);
lambda_omegaQG_RHS_mode = Lambda_z(omega_QG_RHS_mode(:,:,:,:));


figure(2)
plot(lambda_omega_mode,level); hold on;
plot(lambda_omegaQG_mode,level); hold on;
plot(lambda_omegaQG_RHS_mode,level); hold on;
legend('Omega','QG','RHS')
title('Asymmetry and r factor for r=0.7 Mode')
set(gca, 'YDir','reverse')

figure(3)
plot(lambda_omega,level); hold on;
plot(lambda_omegaQG,level); hold on;
plot(lambda_omegaQG_RHS,level); hold on;
legend('Omega','QG','RHS')
title('Asymmetry and r factor for r=0.7 Turbulent')
set(gca, 'YDir','reverse')

end

figure(5)
contour(omega(:,:,12,end)); colorbar
caxis([-2 2])

figure(6)
contour(omega_QG(:,:,12,end)); colorbar
caxis([-2 2])

figure(7); plot(omega(15,:,12,end)); hold on; plot(omega_QG(15,:,12,end))

figure(8)
contour(omega_mode(:,:,12,end)); colorbar

figure(9)
contour(omega_QG_mode(:,:,12,end)); colorbar

figure(10)
plot(omega_mode(15,:,12,end)); hold on;
plot(omega_QG_mode(15,:,12,end))
