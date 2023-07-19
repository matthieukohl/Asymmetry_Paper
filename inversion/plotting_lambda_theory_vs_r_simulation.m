% plot lambda theory vs. turbulent simulation

clear; close all

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
g = 9.80665;  
epsilon = 0.621; % Mass ratio Water Vapour / Air


list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

[lambda_omega,lambda_omega_wzm,lambda_toy_mid,lambda_toy_500,...
    lambda_toy_mid_rhs,lambda_toy_500_rhs, k_av_r,k_mid_r,k_500_r,...
    k_av_r_rhs,k_mid_r_rhs,k_500_r_rhs] = deal(zeros(size(R_factor)));

tic
for tt = 1:10

tt    
    
% Define Reduction Factor

r=R_factor(tt); 

% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

k_av = 0;
k_mid = 0;
k_500 = 0;
k_av_rhs = 0;
k_mid_rhs = 0;
k_500_rhs = 0;
counter = 0;
counter_index = 0;

lambda_place = 0;
lambda_wzm_place = 0;

for i= 9:9
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/ug_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/vg_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');

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

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'ug', start, count);
v_in= ncread(path_v,'vg', start, count);
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
    
% % varying f    
% f = 2 * Omega * sind(lat);
% beta = 1 / R * d_dphi(f, dphi);

% constant f0
lat0 = lat(1)+(lat(end)-lat(1))/2;
f = 2 * Omega * sind(lat0*ones(size(lat)));
beta = 2 * Omega * cosd(lat0*ones(size(lat)))/R;


% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta1 = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

err_rms = zeros(length(100),1);
for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);


sigma_accu(:,:,:,t) = -Ra * T(:,:,:,t)./ (Level .* theta) .* d_dp(theta, level);



% term A (Q vector)

for k = 1 : length(level)
    
    du_dlambda(:, :, k, t) = d_dlambda(u(:, :, k, t), dlambda);
    dv_dlambda(:, :, k, t) = d_dlambda(v(:, :, k, t), dlambda);
    du_dphi(:, :, k, t) = d_dphi(u(:, :, k, t), dphi);
    dv_dphi(:, :, k, t) = d_dphi(v(:, :, k, t), dphi);

    dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
    dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
    Qx(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
        (1 ./ cos(Phi))    .* dv_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
    Qy(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                              dv_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
    div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
    A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
end

clear('temp_x', 'temp_y', 'div_temp');

% term B (planetary vorticity)

temp = d_dp(v(:, :, :, t), level);
for j = 1 : length(lat)
    B(j, :, :, t) = f(j) * beta(j) * temp(j, :, :);
end
clear('temp');

% right hand side

rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t); 


% find troposphere
density = rho(T(:,:,:,t),Level);
lapse = -g*density.*d_dp(T(:,:,:,t),level)*10^3; % K/km

% find troposphere based on 2 K/km criteria

lapse_mean = mean(mean(lapse,1),2);

trop_new = find(lapse_mean>-2,1);

% find wavenumbers

[~,p500] = min(abs(level-50000));
[~,plow] = min(abs(level-92500));
surf = 1;
[~,pmid] = min(abs(level-(level(surf)-(level(surf)-level(trop_new))/2)));
[~,p800] = min(abs(level-80000));
[~,p200] = min(abs(level-20000));

[k_av_new,spec_av_new] = Spectrum1(omega(:,:,p800:p200,t),lat,lon);
[k_mid_new,spec_mid_new] = Spectrum1(omega(:,:,pmid,t),lat,lon);
[k_500_new,spec_500_new] = Spectrum1(omega(:,:,p500,t),lat,lon);

[k_av_rhs_new,spec_av_rhs_new] = Spectrum1(rhs(:,:,p800:p200,t),lat,lon);
[k_mid_rhs_new,spec_mid_rhs_new] = Spectrum1(rhs(:,:,pmid,t),lat,lon);
[k_500_rhs_new,spec_500_rhs_new] = Spectrum1(rhs(:,:,p500,t),lat,lon);



[L_mid_new,L_500_new,L_av_new] = Deformation_Radius(sigma_accu(:,:,:,t),f,level,lat,trop_new);

% k_av_new = k_av_new*L_av_new;
% k_mid_new = k_mid_new*L_mid_new;
% k_500_new = k_500_new*L_500_new;
% 
% k_av_rhs_new = k_av_rhs_new*L_av_new;
% k_mid_rhs_new = k_mid_rhs_new*L_mid_new;
% k_500_rhs_new = k_500_rhs_new*L_500_new;

Ld = 1e6/sqrt(2);

k_av_new = k_av_new*Ld;
k_mid_new = k_mid_new*Ld;
k_500_new = k_500_new*Ld;

k_av_rhs_new = k_av_rhs_new*Ld;
k_mid_rhs_new = k_mid_rhs_new*Ld;
k_500_rhs_new = k_500_rhs_new*Ld;

k_av = k_av + k_av_new;
k_mid = k_mid + k_mid_new;
k_500 = k_500 + k_500_new;

k_av_rhs = k_av_rhs + k_av_rhs_new;
k_mid_rhs = k_mid_rhs + k_mid_rhs_new;
k_500_rhs = k_500_rhs + k_500_rhs_new;


counter = counter + 1;

end

omega_wzm = omega - mean(omega,2);

lambda_place = lambda_place + mean(Lambda_weighted(omega(:,:,1:16,:),lat));

lambda_wzm_place = lambda_wzm_place + mean(Lambda_weighted(omega_wzm(:,:,1:16,:),lat));

counter_index = counter_index + 1;

end

% wavenumbers vs. r

k_av_r(tt) = k_av/counter;
k_mid_r(tt) = k_mid/counter;
k_500_r(tt) = k_500/counter;

k_av_r_rhs(tt) = k_av_rhs/counter;
k_mid_r_rhs(tt) = k_mid_rhs/counter;
k_500_r_rhs(tt) = k_500_rhs/counter;


% lambda

%omega = omega - mean(omega,2);

%lambda_omega(tt) = mean(Lambda_weighted(omega(:,:,1:16,:),lat));

lambda_omega(tt) = lambda_place/counter_index;
lambda_omega_wzm(tt) = lambda_wzm_place/counter_index;

[~,lambda_toy_mid(tt)] = toy_model(R_factor(tt),k_mid_r(tt));

[~,lambda_toy_500(tt)] = toy_model(R_factor(tt),k_500_r(tt));

[~,lambda_toy_mid_rhs(tt)] = toy_model(R_factor(tt),k_mid_r_rhs(tt));

[~,lambda_toy_500_rhs(tt)] = toy_model(R_factor(tt),k_500_r_rhs(tt));

end

% plot nondimensionalized wavenumbers

figure(1)
semilogx(R_factor,k_mid_r,'b','linewidth',1.2); hold on;
semilogx(R_factor,k_500_r,'b--','linewidth',1.2); hold on;
semilogx(R_factor,k_mid_r_rhs,'r','linewidth',1.2); hold on;
semilogx(R_factor,k_500_r_rhs,'r--','linewidth',1.2); hold on;
semilogx(R_factor,k_av_r,'b-.','linewidth',1.2); hold on;
semilogx(R_factor,k_av_r_rhs,'r-.','linewidth',1.2);
xlabel('r');
ylabel('Wavenumber (ND)')
set(gca,'FontSize',14)

figure(2)
semilogx(R_factor,lambda_omega,'k','linewidth',1.2); hold on;
semilogx(R_factor,lambda_omega_wzm,'k--','linewidth',1.2); hold on;
semilogx(R_factor,lambda_toy_mid,'b','linewidth',1.2); hold on;
semilogx(R_factor,lambda_toy_500,'b--','linewidth',1.2); hold on;
semilogx(R_factor,lambda_toy_mid_rhs,'r','linewidth',1.2); hold on;
semilogx(R_factor,lambda_toy_500_rhs,'r--','linewidth',1.2)
xlabel('r');
ylabel('\lambda')
set(gca,'FontSize',14)
set(gca,'box','off')
legend('w','wzm','toy mid spec w','toy 500 spec w','toy mid spec rhs','toy 500 spec rhs','Location','SouthWest');
legend boxoff
%saveas(gcf,'figures_paper/lambda_theory_vs_r_simulation','epsc');

