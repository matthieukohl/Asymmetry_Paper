% plot Skew(p) of the rhs

clear; close all;

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

% path

path_rhs = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/rhs_2017_2017_vs_p.nc';



% Read in Grid-Variables Lat/Lon 

lat = ncread(path_rhs,'latitude');
lon = ncread(path_rhs,'longitude');
level = ncread(path_rhs,'level'); level = double(level);
time = ncread(path_rhs,'time');

dphi = abs((lat(2) - lat(1))) / 180.0 * 3.1415926;
dlambda = abs((lon(2) - lon(1))) / 180.0 * 3.1415926;

lat = flip(lat);
level = flip(level);
level = level*1e2;


[Lon,Lat] = meshgrid(lon,lat);

mask = landmask(Lat,Lon);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

years = 2017:1:2017; % 

month = {'jan','mar','may','jul','sep','nov'};

[skew_rhs,skew_rhs_ocean,skew_rhs_land] = deal(zeros(length(level),1));

counter = 0;

for tt = 1:length(years)
 for m = 1:length(month) 
m
time_range = season_selector(path_rhs, month{m}, years(tt), years(tt));

start = [1, 1, 1, time_range(1)];
count = [length(lon), length(lat),length(level),length(time_range)];

T_in= ncread(path_rhs,'t', start, count);
u_in= ncread(path_rhs,'u', start, count);
v_in= ncread(path_rhs,'v', start, count);

u = rearrange(u_in,lat,lon,level,time_range);
v = rearrange(v_in,lat,lon,level,time_range);
T = rearrange(T_in,lat,lon,level,time_range);


% constant f0
lat0 = lat(1)+(lat(end)-lat(1))/2;
f = 2 * Omega * sind(lat0*ones(size(lat)));
beta = 2 * Omega * cosd(lat0*ones(size(lat)))/R;

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(lat), length(lon), length(level), length(time_range)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(time_range)));
[r,rhs, omega_b] = deal(zeros(length(lat), length(lon), length(level), length(time_range)));
[omega_QG,omega_QG_RHS] = deal(nan(length(lat),length(lon), length(level), length(time_range)));

for t = 1 : length(time_range)
t
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(lat), length(lon), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);

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

% for j = 1 : length(lat)
%     B(j, :, :, t) = f0 * beta * temp(j, :, :);
% end
% clear('temp');

% right hand side

rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t);


end

clear('u','v','T');

[a,ind45] = min(abs(lat-45));
[a,p900] = min(abs(level-900*1e2));
[a,p200] = min(abs(level-200*1e2));

Mask = repmat(mask,[1 1 length(level) length(time_range)]);
rhs_ocean = rhs;
rhs_land = rhs;
rhs_ocean(Mask==1) = NaN;
rhs_land(Mask==0) = NaN;


for kk = 1:length(level)  
skew_rhs(kk) = skew_rhs(kk) + nanmean(Skewness(-rhs(:,:,kk,:)));
skew_rhs_ocean(kk) = skew_rhs_ocean(kk) + nanmean(Skewness(-rhs_ocean(:,:,kk,:)));
skew_rhs_land(kk) = skew_rhs_land(kk) + nanmean(Skewness(-rhs_land(:,:,kk,:)));
end

counter = counter + 1;

end

 
end

skew_rhs = skew_rhs/counter;
skew_rhs_ocean = skew_rhs_ocean/counter;
skew_rhs_land = skew_rhs_land/counter;

figure(1)
plot(skew_rhs,level/100,'k'); hold on;
plot(skew_rhs_ocean,level/100,'b'); hold on;
plot(skew_rhs_land,level/100,'r');
set(gca,'Ydir','reverse');
xlabel('Skew rhs')
ylabel('Pressure (hPa)')
xlim([-8 8])
legend('ocean+land','ocean','land'); legend boxoff

% plot some longitudinal samples at 850hPa;

[a,p850]  = min(abs(level-850));
figure('Renderer', 'painters', 'Position', [10 10 1000 300])

%%for tt= 1:size(rhs,4)
tt = 26;
plot(lon,squeeze(rhs_ocean(ind45,:,p850,tt)),'b'); hold on;
plot(lon,squeeze(rhs_land(ind45,:,p850,tt)),'r')
xlabel('lon')
xlim([lon(1) lon(end)])

skew_ocean = skewness(rhs_ocean(ind45,:,p850,tt));
skew_land = skewness(rhs_land(ind45,:,p850,tt));
legend(['ocean=',num2str(round(skew_ocean,2))],['land=',num2str(round(skew_land,2))]); legend boxoff
title(num2str(tt))
%pause(1)

%hold off

%end

