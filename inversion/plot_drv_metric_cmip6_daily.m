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

for tt = 1:length(year_start)

% Define Path 

% daily

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

lat = lat(lat40:lat60);

dphi = abs(lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = abs(lon(2) - lon(1)) / 180.0 * 3.1415926;

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

%event_latspan = [1:length(lat)]; %83 -111
event_latspan = [lat40:lat60];
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

for t = 1:length(event_timespan)

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), t];
count = [length(event_lonspan), length(event_latspan), length(level_indices),1];

start_sf = [event_lonspan(1), event_latspan(1), 1];
count_sf = [length(event_lonspan), length(event_latspan), length(event_timespan)];

u_in= ncread(path_u,'ua', start, count);
v_in= ncread(path_v,'va', start, count);
omega_in= ncread(path_w,'wap', start, count);
T_in= ncread(path_T,'ta', start, count);


for k = 1:length(level)
u(:,:,k) = u_in(:,:,k)';
v(:,:,k) = v_in(:,:,k)';
T(:,:,k) = T_in(:,:,k)';
omega(:,:,k) = omega_in(:,:,k)';
end

clear('u_in','v_in','omega_in','T_in');


f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);

        
[S,Theta,PV_Ertel,theta_dot,heating_r,...
 PV_rate_r,PV_rate_hoskins,PV_adv_z_full,PV_adv_z_wzm,PV_adv_y,...
 PV_p_full,PV_p_wzm,sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level),1));
[r,density,lapse] = deal(zeros(length(event_latspan), length(event_lonspan), length(level)));


% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :) .* (1e5 ./ Level).^ kappa;
sigma_accu(:, :, :) = -Ra * T(:, :, :) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
T_avg(:) = reshape(mean(mean(T(:, :, :), 1), 2), size(level));
theta_avg (:) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:) =  - Ra * T_avg(:) ./ (level .* theta_avg(:)) .* gradient(theta_avg(:), level);
S(:,:,:) = -T(:,:,:)./(theta).*d_dp(theta,level);
Theta(:,:,:) = theta;

% T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
% theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
% sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r(:,:,:) = r_parameter1(T(:,:,:),Level,level);

density(:,:,:) = rho(T(:,:,:),Level);

lapse(:,:,:) = -g*density(:,:,:).*d_dp(T(:,:,:),level)*10^3; % K/km

for k = 1:length(level)
    zeta(:,:,k) = curl(phi, u(:,:,k), v(:,:,k), dphi, dlambda);
end

PV_Ertel(:,:,:) = Ertel_PV(u(:,:,:),v(:,:,:),T(:,:,:)...
    ,f,lat,lon,level,dphi,dlambda,Phi,phi,kappa,p0,g);

% Calculate heating and diabatic PV-generation

% using the r-method

rr = r; rr(omega>0) = 1;

heating_r(:,:,:) = -(1-rr(:,:,:)).*omega(:,:,:).*S(:,:,:);
heating_r(:,:,:) = heating_r(:,:,:).*(1e5./Level).^kappa;

clear('rr');

F = repmat(f,1,size(omega,2),size(omega,3));

PV_rate_r(:,:,:) = -g*(zeta(:,:,:)+F(:,:,:)).*d_dp(heating_r(:,:,:),level);

PV_rate_r(:,:,:) = PV_rate_r(:,:,:)*10^6*3600; % convert to PVU/h

% using the moisture tendency from condensation (ignoring convection)

% dt_tg_condensation = dt_qg_condensation(:,:,:,t)*L/cp;
% theta_dot(:,:,:,t) = -dt_tg_condensation.*(1e5./Level).^kappa;

% with conditional heating

%theta_dot(:,:,:,t) = heating_r(:,:,:,t);

% without conditional heating

%theta_dot = -(1-r(:,:,:,t)).*omega(:,:,:,t).*S(:,:,:,t).*(1e5./Level).^kappa;

% Calculate the RHS of Eq 74a (Hoskins) but expand to avoid dividing by PV
% values that are close to zero:
% PV_dot = PV^2 x d_dtheta(theta_dot/PV) = PV x dtheta_dot_dp/dtheta_dp
% - theta_dot dPV_dp /dtheta_dp

% for ii = 1:length(lat)
%     for jj = 1:length(lon)
%         
% %PV_rate_hoskins(ii,jj,:,t) = PV_Ertel(ii,jj,:,t).^2.*d_dp(theta_dot(ii,jj,:)./PV_Ertel(ii,jj,:,t),level)./d_dp(theta(ii,jj,:),level);
% PV_rate_hoskins(ii,jj,:,t) = PV_Ertel(ii,jj,:,t).*d_dp(theta_dot(ii,jj,:,t),level)./d_dp(theta(ii,jj,:),level)...;
%      -theta_dot(ii,jj,:,t).*d_dp(PV_Ertel(ii,jj,:,t),level)./d_dp(theta(ii,jj,:),level);
%     end 
% end
% 
% dt = 6*3600;
% 
% % 
% PV_rate_r(:,:,:,t) = PV_Ertel(:,:,:,t).*d_dp(theta_dot(:,:,:,t),level)./d_dp(theta,level);
% 
% PV_rate_r(:,:,:,t) = PV_rate_r(:,:,:,t).*3600;

% calculate the DRV metric

PV_Ertel_prime = PV_Ertel-mean(PV_Ertel,2);

low = 1;
top = length(level);

metric_drv = PV_Ertel_prime.*PV_rate_r;
metric_drv = squeeze(trapz(level(low:top),metric_drv(:,:,low:top,:),3));
metric_drv(metric_drv<0) = 0;

metric_eddy = PV_Ertel_prime.*PV_Ertel_prime;
metric_eddy = squeeze(trapz(level(low:top),metric_eddy(:,:,low:top,:),3));

metric_tend = PV_rate_r.*PV_rate_r;
metric_tend = squeeze(trapz(level(low:top),metric_tend(:,:,low:top,:),3));

metric_drv_new = metric_drv_new + max(max(metric_drv.^2,[],1),[],2);

%metric_drv_vs_surf(tt) = squeeze(mean(metric_drv_vs_t));

metric_drv_norm_new = metric_drv_norm_new + max(max(metric_drv.^2,[],1),[],2)./max(max(metric_eddy.*metric_tend,[],1),[],2);

%metric_vs_surf(tt) = squeeze(mean(metric_vs_t));

metric_drv_norm_x_new = metric_drv_norm_x_new + mean(max(metric_drv.^2,[],2)./max(metric_eddy.*metric_tend,[],2),1);

%metric_x_vs_surf(tt) = squeeze(mean(metric_x_vs_t));


counter = counter + 1;
end

metric_drv_abs(tt) = metric_drv_new/counter;
metric_drv_norm(tt) = metric_drv_norm_new/counter;
metric_drv_norm_x(tt) = metric_drv_norm_x_new/counter;

end







% % plot omega
% cmax = 1;
% nint = 20;
% cint = cmax/nint;
% 
% [a,p500] = min(abs(level-50000));
% 
% figure('Renderer', 'painters', 'Position', [10 10 1200 400])
% 
% contourf(omega(:,:,p500,end-1),-cmax:cint:cmax,'EdgeColor','none'); hold on;
% colorbar
% colormap(redblue(nint));
% caxis([-cmax cmax]);
% title('\rm Vertical Velocity')