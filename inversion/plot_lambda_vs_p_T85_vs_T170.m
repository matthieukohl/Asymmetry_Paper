% compare the w spec from abs T85 to T170

% T85

% plot spectrum of rhs

clear;

list = {'abs1.0_T85','abs4.0_T85'};

lambda_vs_p_T85 = zeros(30,2);

tic
for tt = 1:2
tt    
    
% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

spec_av = 0;
spec_mid = 0;
spec_500 = 0;
k_av = 0;
k_mid = 0;
k_500 = 0;
S_500 = 0;

skew_rhs_500 = 0;
skew_rhs45_500 = 0;
lambda_rhs_500 = 0;

counter = 0;

for i= 1:6 %9:12
i

% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/ug_day',num2str(5 *(28+i)),'0h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/vg_day',num2str(5 *(28+i)),'0h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');


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

event_latspan = [83:111]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [1:200];
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

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        %u(:, :, k, t) = u_in(:, :, k, t)';
        %v(:, :, k, t) = v_in(:, :, k, t)';
        omega(:, :, k, t) = omega_in(:, :, k, t)';
    end
end
    

lambda_vs_p_T85(:,tt) = lambda_vs_p_T85(:,tt) + Lambda_weighted_vs_p(omega,lat);




counter = counter + 1;
end

lambda_vs_p_T85(:,tt) = lambda_vs_p_T85(:,tt)/counter;


end



toc



% T170

list = {'abs1.0_T170','abs4.0_T170'};

lambda_vs_p_T170 = zeros(30,2);

tic
for tt = 1:2

tt    
    
% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

spec_av = 0;
spec_mid = 0;
spec_500 = 0;
k_av = 0;
k_mid = 0;
k_500 = 0;
S_500 = 0;

skew_rhs_500 = 0;
skew_rhs45_500 = 0;
lambda_rhs_500 = 0;

counter = 0;

for i= 1:59 %9:12
i

% Define Path 

path_temp = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/temp_day',num2str(5 *(521+i)),'h00.nc');
path_u = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/u_day',num2str(5 *(521+i)),'h00.nc');
path_v = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/v_day',num2str(5 *(521+i)),'h00.nc');
path_omega = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/omega_day',num2str(5 *(521+i)),'h00.nc');
path_dt_qg_condensation = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/dt_qg_condensation_day',num2str(5 *(521+i)),'h00.nc');



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
%event_latspan = [94:98]; %
event_lonspan = [1:512]; % 1:256
event_timespan = [1:20];
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

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        %u(:, :, k, t) = u_in(:, :, k, t)';
        %v(:, :, k, t) = v_in(:, :, k, t)';
        omega(:, :, k, t) = omega_in(:, :, k, t)';
    end
end
    
lambda_vs_p_T170(:,tt) = lambda_vs_p_T170(:,tt) + Lambda_weighted_vs_p(omega,lat);



counter = counter + 1;

end

lambda_vs_p_T170(:,tt) = lambda_vs_p_T170(:,tt)/counter;

end



% plot lambda vs p

low = 4;

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)

plot(lambda_vs_p_T85(low:end,1),level(low:end)/100,'b'); hold on;
plot(lambda_vs_p_T170(low:end,1),level(low:end)/100,'r'); hold on;
xlabel('$\lambda$')
ylabel('Pressure (hPa)');
set(gca,'Ydir','reverse')
title('\rm abs1.0')
legend('T85','T170'); legend boxoff

subplot(1,2,2)

plot(lambda_vs_p_T85(low:end,2),level(low:end)/100,'b'); hold on;
plot(lambda_vs_p_T170(low:end,2),level(low:end)/100,'r'); hold on;
xlabel('$\lambda$')
ylabel('Pressure (hPa)');
set(gca,'Ydir','reverse');
title('\rm abs4.0')
legend('T85','T170'); legend boxoff

load('figures_paper/trop_height.mat');
strat_level = trop_height;
strat_level = strat_level(14);

nanmean(lambda_vs_p_T85(low:strat_level,2))
nanmean(lambda_vs_p_T170(low:strat_level,2))




