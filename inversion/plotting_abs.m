% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (abs-simulations) iteratively using Dirichlet-BC w=0 on 
% the boundaries.

clear;

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
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

for i=1:16

%i=16;

% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.01/u_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.01/v_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.01/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');

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
event_timespan = [1:100];
level_indices = [5:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'u_plevel', start, count);
v_in= ncread(path_v,'v_plevel', start, count);
omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
    for k = 1 : length(level_indices)
        for t = 1 : length(event_timespan)
            T(:, :, k, t) = T_in(:, :, k, t)';
            u(:, :, k, t) = u_in(:, :, k, t)';
            v(:, :, k, t) = v_in(:, :, k, t)';
            omega(:, :, k, t) = omega_in(:, :, k, t)';
        end
    end
    
    
f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);
            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)
    
% Smooth the T-fields

%deal with points where sigma_accu is not smooth
sigma_smooth_window = 9;    
if sigma_smooth_window ~= 0
    for k = 1 : length(level)
        T(:, :, k, t) = smooth2a(T(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
    end
else
    T(:, :, :, t) = T(:, :, :, t);
end

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

% p_sat - saturation pressure water vapor

p_sat = p_triple*exp(-L/Rc*(1./T(:,:,:,t)-1/T_triple));

% r_sat - saturated mixing ratio

r_sat = epsilon*p_sat./(Level-p_sat); 
 
% Equivalent Potential Temperature (& Gradient)

theta_star = theta.*exp(L/cpc*r_sat./T(:,:,:,t));

dTheta_star_dp = d_dp(theta_star,level);

dTheta_dp = d_dp(theta,level);

% Dry-lapse rate

Gamma_dry = Ra/cp;

% Moist-lapse rate

Gamma_moist = Ra/cp*(1+L/(Ra*T(:,:,:,t)).*r_sat)./(1+(cpc/cp+(L./(Rc*T(:,:,:,t))-1)*L./(cp*T(:,:,:,t))).*r_sat);

% reduction factor r

r = (theta./theta_star).*(Gamma_moist./Gamma_dry).*(dTheta_star_dp./dTheta_dp);

figure(1); contour(r(:,:,12)); colorbar;

r_av(t) = mean(mean(r(:,:,12),1),2);

end

Reduction(i) = mean(r_av)

end




