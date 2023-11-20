clear;

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

omega_QG = [];
omega_QG_RHS = [];
omega = [];
r = [];



for i=1:1
    s = load(strcat('output/abs_mode_simulation/abs6.0_T85/Inversion_smooth5_day.mat'));
    t = load(strcat('output/abs_mode_simulation/abs6.0_T85/Inversion_RHS_smooth5_day.mat'));
    omega_QG=cat(4,omega_QG,s.omega_QG);
    omega_QG_RHS = cat(4,omega_QG_RHS,t.omega_QG);
    omega = cat(4,omega,s.omega);
    r = cat(4,r,t.r);
end

r = r(11:26,:,:,:);

omega_QG = omega_QG(11:26,:,:,:);
omega_QG = omega_QG-mean(omega_QG(:,:,:,:),2);
skewness_omegaQG = Skewness(omega_QG);
lambda_omegaQG = Lambda(omega_QG);

omega_QG_RHS = omega_QG_RHS(11:26,:,:,:);
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
skewness_omegaQG_RHS = Skewness(omega_QG_RHS);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS);

omega = omega(11:26,:,:,:);
omega = omega - mean(omega(:,:,:,:),2);
skewness_omega = Skewness(omega);
lambda_omega = Lambda(omega);


%% Plotting

figure(1)

plot(squeeze(skewness_omega)); hold on
plot(squeeze(skewness_omegaQG)); hold on
plot(squeeze(skewness_omegaQG_RHS)); hold on;
title('Skewness Vertical Velocity (calculated r) abs=4.0')
legend('omega','omegaQG_{total}','omegaQG_{RHS}');
xlabel('time(quarter days)');
ylabel('skewness')

figure(2)

plot(squeeze(lambda_omega)); hold on
plot(squeeze(lambda_omegaQG)); hold on
plot(squeeze(lambda_omegaQG_RHS)); hold on;
title('Lambda Vertical Velocity (calculated r) abs=4.0')
legend('omega','omegaQG_{total}','omegaQG_{RHS}');
xlabel('time(quarter days)');
ylabel('Lambda')

figure(4)
contour(omega(:,:,12,5)); colorbar

figure(5)
contour(omega_QG(:,:,12,5)); colorbar

figure(6)
for tt = 1:15
    
    contour(omega_QG(:,:,12,tt)); colorbar;
    pause(0.5)
end

figure(7)
for tt = 1:15
    
    contour(r(:,:,12,tt)); colorbar;
    pause(0.5)
end





