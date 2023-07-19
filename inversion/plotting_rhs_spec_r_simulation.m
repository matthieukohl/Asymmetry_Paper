% plot spectrum of rhs

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

[spec_av_r,spec_mid_r,spec_500_r] = deal(zeros(10,256));
[spec_500_r_vs_lat] = deal(zeros(10,29,256));
[k_av_r,k_mid_r,k_500_r,S_500_r,lambda_rhs_500_r,skew_rhs_500_r,skew_rhs45_500_r] = deal(zeros(length(R_factor),1));
S_500_r_vs_lat = zeros(10,29);


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

spec_av = 0;
spec_mid = 0;
spec_500 = 0;
spec_500_vs_lat = 0;
k_av = 0;
k_mid = 0;
k_500 = 0;
S_500 = 0;
S_500_vs_lat = zeros(29,1);
skew_rhs_500 = 0;
skew_rhs45_500 = 0;
lambda_rhs_500 = 0;

counter = 0;

for i= 9:9 %9:12
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
%event_latspan = [94:98]; %
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


% [k_av_new,spec_av_new] = Spectrum1(rhs(:,:,plow:trop_new,t),lat,lon);
% [k_mid_new,spec_mid_new] = Spectrum1(rhs(:,:,pmid,t),lat,lon);
% [k_500_new,spec_500_new] = Spectrum1(rhs(:,:,p500,t),lat,lon);

[k_av_new,spec_av_new] = Spectrum2(rhs(:,:,plow:trop_new,t),lat,lon);
[k_mid_new,spec_mid_new] = Spectrum2(rhs(:,:,pmid,t),lat,lon);
[k_500_new,spec_500_new] = Spectrum2(rhs(:,:,p500,t),lat,lon);
[pholder,spec_500_new_vs_lat] = Spectrum2_vs_lat(rhs(:,:,p500,t),lat,lon);



S_500_new = mean(mean(sigma_accu(:,:,p500,t),1),2);
S_500_new_vs_lat = squeeze(mean(sigma_accu(:,:,p500,t),2));

k_av = k_av + k_av_new;
k_mid = k_mid + k_mid_new;
k_500 = k_500  + k_500_new;

S_500 = S_500 + S_500_new;
S_500_vs_lat = S_500_vs_lat + S_500_new_vs_lat;

spec_av = spec_av + spec_av_new;
spec_mid = spec_mid + spec_mid_new;
spec_500 = spec_500  + spec_500_new;
spec_500_vs_lat = spec_500_vs_lat + spec_500_new_vs_lat;

[a,ind] = min(abs(lat-45));

skew_rhs_500 = skew_rhs_500 + Skewness(-rhs(:,:,p500,t));
skew_rhs45_500 = skew_rhs45_500 + skewness(rhs(ind,:,p500,t));
lambda_rhs_500 = lambda_rhs_500 + Lambda_weighted(rhs(:,:,p500,t),lat);


counter = counter + 1;

end

end

spec_av_r(tt,:) = spec_av/counter;
spec_mid_r(tt,:) = spec_mid/counter;
spec_500_r(tt,:) = spec_500/counter;
spec_500_r_vs_lat(tt,:,:) = spec_500_vs_lat/counter;
k_av_r(tt) = k_av/counter;
k_mid_r(tt) = k_mid/counter;
k_500_r(tt) = k_500/counter;
S_500_r(tt) = S_500/counter;
S_500_r_vs_lat(tt,:) = S_500_vs_lat/counter;

skew_rhs_500_r(tt) = skew_rhs_500/counter;
skew_rhs45_500_r(tt) = skew_rhs45_500/counter;
lambda_rhs_500_r(tt) = lambda_rhs_500/counter;

end

toc

R_eff = cosd(lat(12))*R;
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
Ld = sqrt(S_500_r)*800*1e2/1e-4;
Ld = Ld/sqrt(2);
Ld = 0.5*Ld; % because our H is the layer height and not the trospopheric height
% Ld = 1e6/sqrt(2);
% 
% k = k*Ld(1);

figure(1)
semilogx(k(1:N/2)*Ld(1),spec_500_r(1,1:N/2)); hold on;
ylabel('Magnitutde')
xlabel('Wavenumber k')
title('\rm Power spectrum rhs at 500hPa for r=0.01')


k_nd = k.*Ld;
spec = spec_500_r(1,:);
spec = spec./(1+k_nd.^2).^2;

sum(k_nd(1:N/2).*spec(1:N/2))./(sum(spec(1:N/2)))



save('spec_rhs_gcm_r001.mat','k','spec_500_r');

%save('spec_rhs_gcm_vs_lat.mat','lat','lon','k','spec_500_r_vs_lat','S_500_r_vs_lat');



figure(2)
semilogx(R_factor,k_500_r.*Ld); hold on;
semilogx(R_factor,k_500_r.*0.5*1e6/sqrt(2));
xlabel('r')
ylabel('Nondimensional k')
legend('Ld calculated','Ld estimated')
legend boxoff

figure(1)
semilogx(k*Ld(1),spec_500_r(1,:)/max(spec_500_r(1,:))); hold on;
xlabel('k')
legend('r=0.01')
legend boxoff
title('\rm spec rhs 500hPa')

figure(2)
loglog(k*Ld(1),spec_500_r(1,:)/max(spec_500_r(1,:))); hold on;
xlabel('k')
legend('r=0.01')
legend boxoff
title('\rm spec rhs 500hPa')

figure(10)
semilogx(k*Ld(1),spec_500_r(1:3:end,:)); hold on;
xlabel('k')
legend('r=0.01','r=0.1','r=0.4','r=1.0','k^{-3/2}')
legend boxoff

figure(20)
loglog(k*Ld(1),spec_500_r(1:3:end,:)); hold on;
xlabel('k')
legend('r=0.01','r=0.1','r=0.4','r=1.0','k^{-3/2}')
legend boxoff


figure(1)
loglog(k,spec_mid_r(1:3:end,:)); hold on;
loglog(k,5e-32*k.^(-3/2),'k:')
xlabel('k')
legend('r=0.01','r=0.1','r=0.4','r=1.0','k^{-3/2}')
legend boxoff

figure(2)
semilogx(R_factor,lambda_rhs_500_r);
xlabel('r');
ylabel('\lambda');

figure(3)
semilogx(R_factor,skew_rhs_500_r);
xlabel('r');
ylabel('Skew -rhs');
title('\rm r-simulations')

figure(4)
semilogx(R_factor,skew_rhs45_500_r);
xlabel('r');
ylabel('Skew -rhs');
title('\rm r-simulations')

time = 90;

figure('Renderer', 'painters', 'Position', [10 10 900 300]); 
contourf(-omega(:,:,11,time)); colorbar; 
caxis([min(min(omega(:,:,11,time))) -min(min(omega(:,:,11,time)))]); 
colormap redblue

figure('Renderer', 'painters', 'Position', [10 10 900 300]); 
contourf(rhs(:,:,11,time)); colorbar; 
caxis([-max(max(rhs(:,:,11,time))) max(max(rhs(:,:,11,time)))]); 
colormap redblue

