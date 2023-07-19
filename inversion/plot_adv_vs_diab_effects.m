% Solve Hoskins equation 74a to determine the relative contributions from
% vertical adv. vs diabatic effects to the PV tendency in different climates
% Calculate theta_dot by assuming a cosine profile for w and saturated ascent taking the T(p)
% profiles from abs1.0 and abs4.0 respectively.
% d_t PV = dtheta_dot_dp x dp_dtheta x PV + theta_dot x dp_dtheta x dPV_dp
% d_t PV = diab*PV + adv*dPV_dp

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
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air
g = 9.80665;



list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

strat_level_turbulent = [8,9,10,10,11,11,12, 12, 13, 13, 14, 14, 15,16,18];

index = 11*ones(size(surf_temp));

p_low = 5;
p_up = 17;


for tt = 14:14
tt
for i= 6:6
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
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

event_latspan = [90:104]; %83 -111
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

T_in= ncread(path_temp,'temp_plevel', start, count);
omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        T(:, :, k, t) = T_in(:, :, k, t)';
        omega(:, :, k, t) = omega_in(:, :, k, t)';
    end
end
       
f0 = 2 * Omega * sind(45);

[S,Theta,r,theta_dot,dtheta_dot_dp,dp_dtheta,dtheta_dp] = ...
        deal(zeros(size(T)));
 

for t = 1 : length(event_timespan)
    
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);
Theta(:,:,:,t) = theta;

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);

rr = r(:,:,:,t); rr(omega(:,:,:,t)>=0) = 1;

theta_dot(:,:,:,t) = -(1-rr).*omega(:,:,:,t).*S(:,:,:,t);
theta_dot(:,:,:,t) = theta_dot(:,:,:,t).*(1e5./Level).^kappa;

clear('rr');

dtheta_dot_dp(:,:,:,t) = d_dp(theta_dot(:,:,:,t),level);
dp_dtheta(:,:,:,t) = 1./(d_dp(theta,level));
dtheta_dp(:,:,:,t) = d_dp(theta,level);

end

end

diab = dtheta_dot_dp.*dp_dtheta; diab = squeeze(mean(mean(mean(diab(:,:,p_low:p_up,:),1),2),4));

adv = theta_dot.*dp_dtheta; adv = squeeze(mean(mean(mean(adv(:,:,p_low:p_up,:),1),2),4));

q_init = -g*f0*squeeze(mean(mean(mean(dtheta_dp(:,:,p_low:p_up,:),1),2),4));


end

% Integrate the equation

t_end = 2.5*24*3600;

[t,q] = ode45(@(t,q) PV1D(t,q,adv,diab,level(p_low:p_up)),[0 t_end],q_init);

q_prime = q-q_init';

% Make a plot of initial and final PV

lab = list{tt}; lab = lab(1:6);

figure(1)
plot(q(end,:)*1e6,level(p_low:p_up)/100); hold on;
plot(q_init*1e6,level(p_low:p_up)/100)
xlabel('PVU')
ylabel('p(hPa)')
set(gca,'Ydir','reverse')
legend('q(t)','q_{0}'); legend boxoff
ylim([level(p_up)/100 level(p_low)/100])
set(gca,'FontSize',12)
title(num2str(lab))

figure(2)
plot(q_prime(end,:)*1e6,level(p_low:p_up)/100); hold on;
%plot(q_prime(20,:)*1e6,level(p_low:p_up)/100); hold on;
%plot(q_prime(10,:)*1e6,level(p_low:p_up)/100); hold on;
plot(zeros(size(q(end,:))),level(p_low:p_up)/100,'Color','k','Linestyle',':')
xlabel('PVU')
ylabel('p(hPa)')
legend('q_{prime}'); legend boxoff
set(gca,'Ydir','reverse')
title('PV-Perturbations')
ylim([level(p_up)/100 level(p_low)/100])
set(gca,'FontSize',12)
title(num2str(lab))

for ii = 1:size(q,1)
    mag(ii) = rms(q(ii,:))/rms(q(1,:));
    
end

figure(3)
plot(t/(24*3600),mag)
xlabel('t(days)')
ylabel('rms(q(p))')
set(gca,'FontSize',12)
title(num2str(lab))

% Make a vertical PV budget
dt = t(end)-t(end-1);
sigma = 1/dt*log(q(end,:)'./(q(end-1,:)'));

PV_t = sigma.*q(end,:)'; PV_t = PV_t*1e6*3600;
PV_adv = -adv.*squeeze(d_dp(q(end,:)',level(p_low:p_up))); PV_adv = PV_adv*1e6*3600;
PV_diab = diab.*squeeze(q(end,:)'); PV_diab = PV_diab*1e6*3600;
res = PV_t - PV_adv -PV_diab;

figure(4)
plot(PV_t,level(p_low:p_up)/100); hold on;
plot(-PV_adv,level(p_low:p_up)/100); hold on;
title(num2str(lab))
plot(-PV_diab,level(p_low:p_up)/100); hold on;
plot(res,level(p_low:p_up)/100);
xlabel('PVU/h')
ylabel('p(hPa)')
legend('PV_t','PV_{adv}','PV_{diab}','res'); legend boxoff
set(gca,'Ydir','reverse')
ylim([level(p_up)/100 level(p_low)/100])
xlim([-max(PV_t) max(PV_t)])


% % PV-Budget for the perturbations
% 
% PV_t = sigma.*q_prime(end,:)'; PV_t = PV_t*1e6*3600;
% PV_adv = -adv.*squeeze(d_dp(q_prime(end,:)',level(p_low:p_up))); PV_adv = PV_adv*1e6*3600;
% PV_diab = diab.*squeeze(q_prime(end,:)'); PV_diab = PV_diab*1e6*3600;
% res = PV_t - PV_adv -PV_diab;
% 
% figure(5)
% plot(PV_t,level(p_low:p_up)/100); hold on;
% plot(-PV_adv,level(p_low:p_up)/100); hold on;
% title(num2str(lab))
% plot(-PV_diab,level(p_low:p_up)/100); hold on;
% plot(res,level(p_low:p_up)/100);
% xlabel('PVU/h')
% ylabel('p(hPa)')
% legend('PV_t','PV_{adv}','PV_{diab}','res'); legend boxoff
% set(gca,'Ydir','reverse')
% ylim([level(p_up)/100 level(p_low)/100])
% xlim([-max(PV_t) max(PV_t)])





% figure(5)
% plot(diab,level(p_low:p_up)/100); hold on;
% set(gca,'Ydir','reverse')
% legend('diab'); legend boxoff
% ylim([level(p_up)/100 level(p_low)/100])
% 
% figure(6)
% plot(adv,level(p_low:p_up)/100); hold on;
% set(gca,'Ydir','reverse')
% legend('adv'); legend boxoff
% ylim([level(p_up)/100 level(p_low)/100])




