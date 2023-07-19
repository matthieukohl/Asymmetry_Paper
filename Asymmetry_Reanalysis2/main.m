% Calculate the Asymmetry Parameter from ERAI-Reanalysis and Compare to
% Theoretical Prediction

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
p0 = 101325; %US-Standard Atmosphere
L = 2493*10^3; % Latent Heat Water Vapour
g = 9.80665;  
epsilon = 0.621; % Mass ratio Water Vapour / Air
window_x = 1;
window_y = 1;
km = 1000;
day = 24*3600;
dt = 3600; % 1-hour ERA Resolution
%drag = 2;

% choose Hemisphere

hemisphere = 'north';

tf_hemisphere = strcmp(hemisphere,'north');

% choose Data-Set

data = 'erai';

tf_data = strcmp(data,'erai');

% load the data

if tf_data==1

path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai_monthly_mean_2001_2011.nc');

else
    
path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5_monthly_mean_2001_2011_unpack.nc');
   
end

% Read in Grid-Variables Lat/Lon 

lat_series_original_single = ncread(path, 'latitude');
lon_series_original_single = ncread(path, 'longitude');
lat_series_original = double(lat_series_original_single);
lon_series_original  = double(lon_series_original_single);
latmax = 90; latmin = -90; lonmax = 360; lonmin = 0;
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path, 'level')); %flipud(double(ncread(path_temp, 'level')));
lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices); %flipud(lat_series_original(lat_indices));
dphi = abs(lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = abs(lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

if tf_hemisphere ==1

% % Northern Hemisphere: 30-70 degrees

[~,lat30] = min(abs(lat_series_original-30));
[~,lat70] = min(abs(lat_series_original-70));
event_latspan = [lat70:lat30]; 
event_lonspan = [1:length(lon_series_original)]; 

else

% Southern Hemisphere: 30-70 degrees

[~,latm30] = min(abs(lat_series_original-(-30)));
[~,latm70] = min(abs(lat_series_original-(-70)));
event_latspan = [latm30:latm70]; 
event_lonspan = [1:length(lon_series_original)]; 

end

level_indices = [1:37]; 
time_range = 1:1:12*11;

level = 100*plevels(level_indices);
lat = lat_series(event_latspan);
lon = lon_series(event_lonspan);


% flip lat (to go from low to high lat)
lat = flip(lat);
[Lon,Lat] = meshgrid(lon,lat);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);
% flip level (to go from surface to top)
level = flip(level);

% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

% Read in Variables (Lat/Lon need to be flipped)

[T_in,omega_in] = deal(zeros(length(event_lonspan), length(event_latspan), length(level_indices), length(time_range)));


for tt = 1:length(time_range)
    
start = [event_lonspan(1), event_latspan(1), level_indices(1), time_range(tt)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), 1];

start_sf = [event_lonspan(1), event_latspan(1), time_range(tt)];
count_sf = [length(event_lonspan), length(event_latspan), 1];

T_in(:,:,:,tt) = ncread(path,'t', start, count);
omega_in(:,:,:,tt) = ncread(path,'w', start, count);

end

omega = rearrange(omega_in,lat,lon,level,time_range);
T = rearrange(T_in,lat,lon,level,time_range);

[lambda_era,lambda_theory,Wavenumber,r_mid,L_mid,T_surf] = deal(zeros(length(size(omega,4)),1));

[sigma_accu,S,r,density,T_virtual,lapse,Theta] = deal(zeros(size(omega)));

for t = 1:length(time_range)

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (p0 ./ Level).^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);
density(:,:,:,t) = rho(T(:,:,:,t),Level);
T_virtual(:,:,:,t) = temp_virtual(T(:,:,:,t),Level);
lapse(:,:,:,t) = -g*density(:,:,:,t).*d_dp(T(:,:,:,t),level)*10^3; % K/km

Theta(:,:,:,t) = theta;

% find troposphere based on 2 K/km criteria

lapse_mean = mean(mean(lapse(:,:,:,t),1),2);

trop = find(lapse_mean>-2,1);

% find the dominant Wavenumber

[Wavenumber(t)] = Spectrum1(omega(:,:,1:trop,t),lat,lon);

% find the Non-Dimensionalization to use

[L_mid(t)] = Deformation_Radius(sigma_accu(:,:,:,t),f,level,lat,trop);

% find the r-value to use

[r_mid(t)] = reduction_factor(r(:,:,:,t),level,trop);

% find the theoretical prediction for the asymmetry

lambda_theory(t) = toy_model(r_mid(t),Wavenumber(t),L_mid(t));


% find lambda from reanalysis


%lambda_era(t) = Lambda_weighted(-omega(:,:,1:trop,t),lat);

lambda_era(t) = Lambda(-omega(:,:,1:trop,t));

% find T_surf

T_surf(t) = mean(mean(T(:,:,1,t),1),2);


end


% Average over the months between years 2001:2011

[lambda_theory_av,lambda_era_av,r_mid_av,Wavenumber_av,L_mid_av,T_surf_av] = deal(zeros(12,1));

for ii = 1:12
    
indx_month = ii:12:length(time_range);
    
lambda_theory_av(ii) = mean(lambda_theory(indx_month));

lambda_era_av(ii) = mean(lambda_era(indx_month));

r_mid_av(ii) = mean(r_mid(indx_month));

L_mid_av(ii) = mean(L_mid(indx_month));

Wavenumber_av(ii) = mean(Wavenumber(indx_month));

T_surf_av(ii) = mean(T_surf(indx_month));

end
    
% Make Plots of Seasonal Cycle
figure(1)
plot(lambda_era_av,'Linewidth',1.5); hold on;
plot(lambda_theory_av,'linewidth',1.5)
legend('ERA','Theory'); legend boxoff
ylabel('\lambda')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)


figure(2)
plot(Wavenumber_av.*L_mid_av,'Linewidth',1.5)
ylabel('k')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)


figure(3)
plot(r_mid_av,'Linewidth',1.5)
ylabel('r')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)


figure(4)
plot(L_mid_av/1000,'Linewidth',1.5)
ylabel('L_D(km)')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)

figure(5)
plot(T_surf_av,'Linewidth',1.5)
ylabel('T(K)')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)





















































