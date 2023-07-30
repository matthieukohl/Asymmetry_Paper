% plot lambda from T170 simulations

clear;
close all;

% Define Constants

Omega = 7.2921e-5; % Rotation Rate
R = 6371000.0; % Planetary Radius 
Ra = 287.04; % Gas-constant Air
Rc = 461.92; % Gas-constant Water Vapour
cp = 1005.7; % Heat Capacity Air
cpc = 1850; % Heat Capacity Water Vapour
kappa = Ra/cp;
T_triple = 273.16; % Triple Point Temperature
p_triple = 611; % Triple Point Pressure
p0 = 101325; %US-Standard Atmosphere
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air
g = 9.80665;     % gravitational acceleration [m/s^2]
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% Plot DRV storms in the macroturbulent simulation

list = {'abs4.0_T170'};
%list = {'abs4.0_T85'};

for tt = 1:1
lambda_omega = [];
counter = 0;
 
for i= 1:59
   
i
% Define Path 

path_temp = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/temp_day',num2str(5 *(521+i)),'h00.nc');
path_u = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/u_day',num2str(5 *(521+i)),'h00.nc');
path_v = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/v_day',num2str(5 *(521+i)),'h00.nc');
path_omega = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/omega_day',num2str(5 *(521+i)),'h00.nc');
path_dt_qg_condensation = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/dt_qg_condensation_day',num2str(5 *(521+i)),'h00.nc');


% path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
% path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/u_day',num2str(5 *(28+i)),'0h00.nc');
% path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/v_day',num2str(5 *(28+i)),'0h00.nc');
% path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');
% path_dt_qg_condensation = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/dt_qg_condensation_day',num2str(5 *(28+i)),'0h00.nc');


% Read in Grid-Variables Lat/Lon 

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
latmax = 90; latmin = -90; lonmax = 360; lonmin = 0;
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));


lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

event_latspan = [186:215]; %83 -111
event_lonspan = [1:512]; % 1:256
%event_latspan = [83:113]; %83 -111
%event_lonspan = [1:256]; % 1:256
event_timespan = [1:20]; %Paul's example
level_indices = [1:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];


omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor,dt_qg_condensation] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        omega(:, :, k, t) = omega_in(:, :, k, t)';
        
    end
end


lambda_pholder = Lambda_vs_lat(omega(:,:,1:17,:));     


lambda_omega = cat(2,lambda_omega,lambda_pholder);



counter = counter + 1;
end

end

lambda_mean = mean(lambda_omega,2);
lambda_max = max(lambda_omega,[],2);
lambda_min =  min(lambda_omega,[],2);

figure(2)
subplot(1,2,2)
histogram(mean(lambda_omega,1),[0.5:0.01:0.9],'normalization','probability');

figure(1)
plot(lat,lambda_mean,'k'); hold on;
plot(lat,lambda_max,'r'); hold on;
plot(lat,lambda_min,'b'); hold on;
xlabel('lat');
ylabel('\lambda');
ylim([0.5 1])




list = {'abs4.0_T85'};

for tt = 1:1
lambda_omega = [];
counter = 0;
 
for i= 1:6
   
i
% Define Path 

% path_temp = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/temp_day',num2str(5 *(521+i)),'h00.nc');
% path_u = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/u_day',num2str(5 *(521+i)),'h00.nc');
% path_v = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/v_day',num2str(5 *(521+i)),'h00.nc');
% path_omega = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/omega_day',num2str(5 *(521+i)),'h00.nc');
% path_dt_qg_condensation = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/dt_qg_condensation_day',num2str(5 *(521+i)),'h00.nc');
% 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/u_day',num2str(5 *(28+i)),'0h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/v_day',num2str(5 *(28+i)),'0h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');
path_dt_qg_condensation = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/dt_qg_condensation_day',num2str(5 *(28+i)),'0h00.nc');


% Read in Grid-Variables Lat/Lon 

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
latmax = 90; latmin = -90; lonmax = 360; lonmin = 0;
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));


lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

event_latspan = [93:108]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [1:200]; %Paul's example
level_indices = [1:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];


omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor,dt_qg_condensation] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        omega(:, :, k, t) = omega_in(:, :, k, t)';
        
    end
end

lambda_pholder = Lambda_vs_lat(omega(:,:,1:17,:));

lambda_omega = cat(2,lambda_omega,lambda_pholder);

counter = counter + 1;
end

end

lambda_mean = mean(lambda_omega,2);
lambda_max = max(lambda_omega,[],2);
lambda_min = min(lambda_omega,[],2);

figure(1)
plot(lat,lambda_mean,'k:'); hold on;
plot(lat,lambda_max,'r:'); hold on;
plot(lat,lambda_min,'b:');
ylim([0.4 1])
xlim([lat(1) lat(end)]);
legend('mean T170','max T170','min T170','mean T85','max T85','min T85');
legend boxoff


figure(2)
subplot(1,2,1)
histogram(mean(lambda_omega,1),[0.5:0.01:0.9],'normalization','probability');




