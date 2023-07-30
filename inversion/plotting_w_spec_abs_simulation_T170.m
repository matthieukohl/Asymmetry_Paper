% plot spectrum of rhs

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

[spec_av_r,spec_mid_r,spec_500_r] = deal(zeros(10,512));
[S_500_r,k_av_r,k_mid_r,k_500_r,lambda_rhs_500_r,skew_rhs_500_r,skew_rhs45_500_r] = deal(zeros(length(R_factor),1));

list = {'abs4.0_T170'};

tic
for tt = 1:1

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

T_in= ncread(path_temp,'temp_plevel', start, count);
%u_in= ncread(path_u,'u_plevel', start, count);
%v_in= ncread(path_v,'v_plevel', start, count);
omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        T(:, :, k, t) = T_in(:, :, k, t)';
        %u(:, :, k, t) = u_in(:, :, k, t)';
        %v(:, :, k, t) = v_in(:, :, k, t)';
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


% find various levels
trop_new = 14;
[~,p500] = min(abs(level-50000));
[~,plow] = min(abs(level-92500));
surf = 1;
[~,pmid] = min(abs(level-(level(surf)-(level(surf)-level(trop_new))/2)));

% find the dominant Wavenumber

%R_eff = cosd(lat(12))*R;
N = length(lon);
%k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);
%k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];

Lnd = 1e6/sqrt(2);
k = k*Lnd;

% [k_av_new,spec_av_new] = Spectrum1(rhs(:,:,plow:trop_new,t),lat,lon);
% [k_mid_new,spec_mid_new] = Spectrum1(rhs(:,:,pmid,t),lat,lon);
% [k_500_new,spec_500_new] = Spectrum1(rhs(:,:,p500,t),lat,lon);

[k_av_new,spec_av_new] = Spectrum2(omega(:,:,plow:trop_new,t),lat,lon);
[k_mid_new,spec_mid_new] = Spectrum2(omega(:,:,pmid,t),lat,lon);
[k_500_new,spec_500_new] = Spectrum2(omega(:,:,p500,t),lat,lon);

S_500_new = mean(mean(sigma_accu(:,:,p500,t),1),2);


k_av = k_av + k_av_new;
k_mid = k_mid + k_mid_new;
k_500 = k_500  + k_500_new;

S_500 = S_500 + S_500_new;

spec_av = spec_av + spec_av_new;
spec_mid = spec_mid + spec_mid_new;
spec_500 = spec_500  + spec_500_new;

[a,ind] = min(abs(lat-45));

counter = counter + 1;

end

end

spec_av_r(tt,:) = spec_av/counter;
spec_mid_r(tt,:) = spec_mid/counter;
spec_500_r(tt,:) = spec_500/counter;
k_av_r(tt) = k_av/counter;
k_mid_r(tt) = k_mid/counter;
k_500_r(tt) = k_500/counter;

S_500_r(tt) = S_500/counter;

skew_rhs_500_r(tt) = skew_rhs_500/counter;
skew_rhs45_500_r(tt) = skew_rhs45_500/counter;
lambda_rhs_500_r(tt) = lambda_rhs_500/counter;

end

toc

R_eff = cosd(lat(12))*R;
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
dp = 800*1e2;
Ld = sqrt(S_500_r)*dp/1e-4;
Ld = Ld/sqrt(2); % because of our choice of nondimensionalization
Ld = 0.5*Ld; % because H is the layer and not the tropospheric height

k = k*Ld(1);

figure(1)
semilogx(k(1:N/2),spec_500_r(1,1:N/2)); hold on;
ylabel('Magnitutde')
xlabel('Wavenumber k')
title('Power spectrum w at 500hPa for r=0.01')


figure(2)
semilogx(R_factor,k_500_r.*Ld); hold on;
semilogx(R_factor,k_500_r.*0.5*1e6/sqrt(2));
xlabel('r')
ylabel('Nondimensional k')
legend('Ld calculated','Ld fixed')
legend boxoff

figure('Renderer', 'painters', 'Position', [10 10 900 400])

subplot(1,2,1)
semilogx(k(1:N/2),spec_500_r(1,1:N/2)/max(spec_500_r(1,1:N/2))); hold on;
ylabel('Magnitutde')
xlabel('Wavenumber k')
title('Power spectrum w at 500hPa for abs4.0 T170')

subplot(1,2,2)
semilogx(k(1:N/2),(k(1:N/2).^2+1).*spec_500_r(1,1:N/2)/(max((k(1:N/2).^2+1).*spec_500_r(1,1:N/2)))); hold on;
ylabel('Magnitutde')
xlabel('Wavenumber k')
title('$(k^2+1)$ Power spectrum w at 500hPa for abs4.0 T170')


