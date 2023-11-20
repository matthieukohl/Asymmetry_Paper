% calculate the r-parameter vs. global mean surface temperature

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
p0 = 101325; %US-Standard Atmosphere
g = 9.80665;     % gravitational acceleration [m/s^2]
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

load('/net/halo/disk28/disk7/mkohl/inversion/figures_paper/trop_height.mat');
strat_level = trop_height;

[r_300,r_500,r_700,r_mid,r_av,r_av_no_bl] = deal(zeros(15,1));

tic
for tt = 1:15
tt

r_300_new = 0;
r_500_new = 0;
r_700_new = 0;
r_mid_new = 0;
r_av_new = 0;
r_av_wbl_new = 0;

counter = 0;


for i= 1:3
%i
% Define Path 

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

event_latspan = [83:113];
event_lonspan = [1:256]; 
event_timespan = [1:200]; 
level_indices = [1:30]; 
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

[T] = deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        T(:, :, k, t) = T_in(:, :, k, t)';
    end
end
       
f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);

        
[S,Theta,PV_Ertel,heating_r,...
 PV_rate_r,PV_rate_hoskins,PV_adv_z_full,PV_adv_z_wzm,PV_adv_y,...
 PV_p_full,PV_p_wzm,sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,density,lapse] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);
Theta(:,:,:,t) = theta;

% T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
% theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
% sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);



% calculate the r-parameter at different levels

[a,p700] = min(abs(level-700*1e2));
[a,p500] = min(abs(level-500*1e2));
[a,p300] = min(abs(level-300*1e2));

r_300_new = r_300_new + nanmean(nanmean(r(:,:,p300,t)));
r_500_new = r_500_new + nanmean(nanmean(r(:,:,p500,t)));
r_700_new = r_700_new + nanmean(nanmean(r(:,:,p700,t)));

r_mid_new = r_mid_new + nanmean(nanmean(r(:,:,round(strat_level(tt)/2),t)));
r_av_new = r_av_new + nanmean(nanmean(nanmean(r(:,:,1:strat_level(tt),t))));
r_av_wbl_new = r_av_wbl_new + nanmean(nanmean(nanmean(r(:,:,4:strat_level(tt),t))));


counter = counter + 1;
end


end

r_300(tt) = r_300_new/counter;
r_500(tt) = r_500_new/counter;
r_700(tt) = r_700_new/counter;

r_mid(tt) = r_mid_new/counter;
r_av(tt) = r_av_new/counter;
r_av_no_bl(tt) = r_av_wbl_new/counter;

end

figure

plot(surf_temp,r_300,'r-o','Linewidth',1.2); hold on;
plot(surf_temp,r_500,'b-o','Linewidth',1.2); hold on;
plot(surf_temp,r_700,'g-o','Linewidth',1.2); hold on;
plot(surf_temp,r_mid,'m-o','Linewidth',1.2); hold on;
plot(surf_temp,r_av,'k-o','Linewidth',1.2);
xlabel('Global Mean Surface Temperature (K)')
ylabel('Reduction Factor r');
set(gca,'FontSize',12)
legend('r_300','r_500','r_700','r_mid','r_av'); legend boxoff;


%save('/net/halo/disk28/disk7/mkohl/inversion/figures_paper/r_vs_global_temp.mat','r_300','r_500','r_700','r_mid','r_av');



% make the figure for the paper

close all; clear;

load('/net/halo/disk28/disk7/mkohl/inversion/figures_paper/r_vs_global_temp.mat','r_500','r_av');

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

figure

plot(surf_temp,r_500,'b','Linewidth',1.2); hold on;
plot(surf_temp,r_av,'k','Linewidth',1.2);
xlabel('Global Mean Surface Temperature (K)')
ylabel('Reduction Factor r');
ylim([0 1])
set(gca,'FontSize',12)
set(gca,'box','off')
legend('500hPa','Tropospheric Average'); legend boxoff;
%saveas(gcf,'/net/halo/disk28/disk7/mkohl/inversion/figures_paper/r_parameter_gcm','epsc')

saveas(gcf,'/net/graupel/archive1/users/mkohl/figures_paper_asymmetry/r_parameter_gcm','epsc');






















