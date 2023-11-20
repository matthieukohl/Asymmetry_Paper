% Plot the DRV metric for CMIP6 runs
% 1%/year C02 increase outputted daily

close all; clear;

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
g = 9.80665;     % gravitational acceleration [m/s^2]
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% generate list of all the files

year_start = 1850:5:2010;
year_end = 1854:5:2014;

[metric_drv_abs,metric_drv_norm,metric_drv_norm_x] = deal(zeros(length(year_start),1));

tic

for ii = 1:2
    
if ii ==1

tt = 1;

else
    
tt = 33;
end

path_u = ['/net/graupel/archive1/users/mkohl/cmip6_daily/ua_day_MPI-ESM1-2-HR_1pctCO2_r1i1p1f1_gn_',num2str(year_start(tt)),'0101-',num2str(year_end(tt)),'1231.nc'];

path_v = ['/net/graupel/archive1/users/mkohl/cmip6_daily/va_day_MPI-ESM1-2-HR_1pctCO2_r1i1p1f1_gn_',num2str(year_start(tt)),'0101-',num2str(year_end(tt)),'1231.nc'];

path_w = ['/net/graupel/archive1/users/mkohl/cmip6_daily/wap_day_MPI-ESM1-2-HR_1pctCO2_r1i1p1f1_gn_',num2str(year_start(tt)),'0101-',num2str(year_end(tt)),'1231.nc'];

path_T = ['/net/graupel/archive1/users/mkohl/cmip6_daily/ta_day_MPI-ESM1-2-HR_1pctCO2_r1i1p1f1_gn_',num2str(year_start(tt)),'0101-',num2str(year_end(tt)),'1231.nc'];


% Read in Grid-Variables Lat/Lon 

lat = ncread(path_w,'lat'); %lat = double(lat); lat = flip(lat);
lon = ncread(path_w,'lon');%lon = double(lon);
level = ncread(path_w,'plev'); %level = double(level)*100; level = flip(level);
time = ncread(path_w,'time');

% Specify Inversion Domain

[a,lat40] = min(abs(lat-40));
[a,lat60] = min(abs(lat-60));

%lat = lat(lat40:lat60);

dphi = abs(lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = abs(lon(2) - lon(1)) / 180.0 * 3.1415926;

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

event_latspan = [1:length(lat)]; %83 -111
%event_latspan = [lat40:lat60];
event_lonspan = [1:length(lon)]; % 1:256
event_timespan = [1:length(time)];
level_indices = [1:length(level)];




% pre-allocation

[u,v,T,omega,A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level)));
   
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

err_rms = zeros(length(event_timespan),1);

metric_drv_new = 0;
metric_drv_norm_new = 0;
metric_drv_norm_x_new = 0;
counter = 0;

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), 1];
count = [length(event_lonspan), length(event_latspan), length(level_indices),length(time)];

start_sf = [event_lonspan(1), event_latspan(1), 1];
count_sf = [length(event_lonspan), length(event_latspan), length(event_timespan)];


T_in= ncread(path_T,'ta', start, count);

for t= 1:length(time)
for k = 1:length(level)
T(:,:,k,t) = T_in(:,:,k,t)';
end

end

clear('T_in');

if ii ==1
    
T_av_old = mean(T(:,:,1,:),4);
else
T_av_new = mean(T(:,:,1,:),4);
end

end

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

[Lon,Lat] = meshgrid(lon,lat);

lon(lon>180) = lon(lon>180)-360;

[Lon_mask,Lat_mask] = meshgrid(lon,lat);

T_av = T_av_new;

mask = landmask(Lat_mask,Lon_mask);
contourf(Lon,Lat,T_av,'EdgeColor','none'); hold on;
contour(Lon,Lat,mask,'k');
colorbar
colormap(flipud(autumn))
title('\rm Surface Air Temperature (K) ')
caxis([min(T_av(:)) max(T_av(:))])

figure
plot(squeeze(mean(T_av_old,2)),lat,'b'); hold on;
plot(squeeze(mean(T_av_new,2)),lat,'r'); hold on;
 xlabel('Temperature (K)')
ylabel('Latitude')
set(gca,'FontSize',12);
ylim([lat(1) lat(end)])
legend('1850','2010'); legend boxoff;

T_av = T_av_new; 
T_av = T_av - 270;
cmax = max(abs(T_av(:)));
nint = 20;
cint = cmax/nint;

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

contourf(Lon,Lat,T_av,-cmax:cint:cmax,'EdgeColor','none'); hold on;
contour(Lon,Lat,mask,'k');
colorbar
colormap(redblue(nint))
title('\rm Surface Air Temperature 2010 (C)')
caxis([-cmax cmax])
xlabel('Longitude')
ylabel('Latitude')



%T_av = T_av - 270;
T_av = T_av_new-T_av_old;
cmax = max(abs(T_av(:)));
nint = 20;
cint = cmax/nint;

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

contourf(Lon,Lat,T_av,-cmax:cint:cmax,'EdgeColor','none'); hold on;
contour(Lon,Lat,mask,'k');
colorbar
colormap(redblue(nint))
title('\rm Surface Air Temperature Anomaly')
caxis([-cmax cmax])
xlabel('Longitude')
ylabel('Latitude')

figure

plot(squeeze(mean(T_av,2)),lat);
ylabel('Latitude')
xlabel('Temp Anomaly (K)')
set(gca,'FontSize',12);
ylim([lat(1) lat(end)])

temps = [mean(mean(T_av_old)),mean(mean(T_av_new))];

bar(temps);
