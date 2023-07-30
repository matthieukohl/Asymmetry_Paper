% Calculate the Asymmetry Parameter from ERA5-Reanalysis and Compare to
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
nmonths = 12;

% choose Hemisphere

hemisphere = 'south';

tf_hemisphere = strcmp(hemisphere,'north');

% choose Data

data = 'erai';
tf_data = strcmp(data,'erai');

% generate list of years

years = 2010:1:2014;
nyears = length(years);
% initiate files

[lambda_era,lambda_era_mid,lambda_theory,k,r_mid,L_mid,T_surf] = ...
    deal(zeros(nyears*nmonths,1));



for kk = 1:length(years)
  years(kk) 
 
for ll = 1:nmonths
ll
tic

% Define Path

% if tf_data==1
% 
% path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai_NH_',num2str(ll),'_',num2str(years(kk)),'.nc');
%    
% else
%     
% path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5_NH_',num2str(ll),'_',num2str(years(kk)),'.nc');
%    
% end 

if tf_data==1
    
    if tf_hemisphere ==1

    path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai/erai_NH_',num2str(ll),'_',num2str(years(kk)),'.nc');
   
    else
      
     path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai/erai_SH_',num2str(ll),'_',num2str(years(kk)),'.nc');
      
    end   
    
else
    if tf_hemisphere ==1
    
    path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_',num2str(ll),'_',num2str(years(kk)),'.nc');
   
    else 
        
    path = strcat('/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ll),'_',num2str(years(kk)),'.nc');
     
    end
    
end 


% Read in Grid-Variables Lat/Lon 

lat_series_original_single = ncread(path, 'latitude');
lon_series_original_single = ncread(path, 'longitude');
lat_series_original = double(lat_series_original_single);
lon_series_original  = double(lon_series_original_single);
latmax = 90; latmin = -90; lonmax = 180; lonmin = -180;
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path, 'level')); %flipud(double(ncread(path_temp, 'level')));
lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices); %flipud(lat_series_original(lat_indices));
dphi = abs(lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = abs(lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

if tf_hemisphere ==1

% % Northern Hemisphere: 30-70 degrees

[~,lat30] = min(abs(lat_series_original-40));
[~,lat70] = min(abs(lat_series_original-60));
event_latspan = [lat70:lat30]; 
event_lonspan = [1:length(lon_series_original)]; 

else

% Southern Hemisphere: 30-70 degrees

[~,latm30] = min(abs(lat_series_original-(-40)));
[~,latm70] = min(abs(lat_series_original-(-60)));
event_latspan = [latm30:latm70]; 
event_lonspan = [1:length(lon_series_original)]; 

end


level_indices = [1:37];

% Define length of time indices taking into account length of month and
% leap year
month31 = [1,3,5,7,8,10,12];
if ismember(ll,month31)
    time_range = 1:1:31*4;
    elseif ll==2
        if mod(years(kk),4)==0
        time_range = 1:1:29*4;
        else
        time_range = 1:1:28*4;
        end
    else
    time_range = 1:1:30*4;
end
   

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

%[T_in,omega_in] = deal(zeros(length(event_lonspan), length(event_latspan), length(level_indices), length(time_range)));
[T,omega] = deal(zeros(length(event_latspan),length(event_lonspan),length(level_indices), length(time_range)));

for tt = 1:length(time_range)
    
start = [event_lonspan(1), event_latspan(1), level_indices(1), time_range(tt)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), 1];

start_sf = [event_lonspan(1), event_latspan(1), time_range(tt)];
count_sf = [length(event_lonspan), length(event_latspan), 1];

T_in = ncread(path,'t', start, count);
omega_in = ncread(path,'w', start, count);

T(:,:,:,tt) = rearrange_era(T_in,lat,lon,level);
omega(:,:,:,tt) = rearrange_era(omega_in,lat,lon,level);
clear('T_in','omega_in')

end

% remove the zonal mean

omega = omega-mean(omega,2);

% omega = rearrange(omega_in,lat,lon,level,time_range);
% T = rearrange(T_in,lat,lon,level,time_range);

%clear('omega_in','T_in')

%[sigma_accu,S,r,density,T_virtual,lapse,Theta] = deal(zeros(size(omega)));
[k_daily,L_mid_daily,r_mid_daily,T_surf_daily,lambda_theory_daily,...
    lambda_era_daily,lambda_era_daily_mid] = deal(zeros(length(time_range),1));
W_power_daily = zeros(size(omega,2),size(omega,4));

for t = 1:length(time_range)

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (p0 ./ Level).^ kappa;
sigma_accu = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

r = r_parameter1(T(:,:,:,t),Level,level);
density = rho(T(:,:,:,t),Level);
lapse = -g*density.*d_dp(T(:,:,:,t),level)*10^3; % K/km

% find troposphere based on 2 K/km criteria

lapse_mean = mean(mean(lapse,1),2);

trop = find(lapse_mean>-2,1);


clear('lapse','theta','density','Level') 

% find the dominant Wavenumber
% 'time spectrum'
% tic
[k_daily(t),spec] = Spectrum1(omega(:,:,1:trop,t),lat,lon);
% toc

W_power_daily(:,t) = spec;

% find the Non-Dimensionalization to use
% 'time ND'
% tic
[L_mid_daily(t)] = Deformation_Radius(sigma_accu,f,level,lat,trop);
% toc
clear('sigma_accu')

% find the r-value to use
% 'time r'
% tic
[r_mid_daily(t)] = reduction_factor(r,level,trop);
% toc
clear('r')

% find the theoretical prediction for the asymmetry
% 'time toy model'
% tic
lambda_theory_daily(t) = toy_model(r_mid_daily(t),k_daily(t),L_mid_daily(t));
% toc
% find lambda from reanalysis
% 'time asymmetry'
% tic

lambda_era_daily(t) = Lambda(omega(:,:,1:trop,t));

surf = 1;
[~,I] = min(abs(level-(level(surf)-(level(surf)-level(trop))/2)));

lambda_era_daily_mid(t) = Lambda(omega(:,:,I,t));

    
% toc
%lambda_era((kk-1)*nmonths+t) = Lambda_weighted(omega(:,:,1:trop,t),lat);

% find T_surf

T_surf_daily(t) = mean(mean(T(:,:,1,t),1),2);

end


clear('omega')
clear('T')

k((kk-1)*nmonths+ll) = mean(k_daily);
r_mid((kk-1)*nmonths+ll) = mean(r_mid_daily);
L_mid((kk-1)*nmonths+ll) = mean(L_mid_daily);
W_power(:,(kk-1)*nmonths+ll) = mean(W_power_daily,2);
lambda_theory((kk-1)*nmonths+ll) = mean(lambda_theory_daily);
lambda_era((kk-1)*nmonths+ll) = mean(lambda_era_daily);
lambda_era_mid((kk-1)*nmonths+ll) = mean(lambda_era_daily_mid);
T_surf((kk-1)*nmonths+ll) = mean(T_surf_daily);

toc
end

end
% Average over the months between years 

[lambda_theory_av,lambda_era_av,lambda_era_mid_av,r_mid_av,k_av,L_mid_av,T_surf_av] = deal(zeros(nmonths,1));

for ii = 1:nmonths
    
indx_month = ii:12:nyears*nmonths;
    
lambda_theory_av(ii) = mean(lambda_theory(indx_month));

lambda_era_av(ii) = mean(lambda_era(indx_month));

lambda_era_mid_av(ii) = mean(lambda_era_mid(indx_month));

r_mid_av(ii) = mean(r_mid(indx_month));

L_mid_av(ii) = mean(L_mid(indx_month));

k_av(ii) = mean(k(indx_month));

T_surf_av(ii) = mean(T_surf(indx_month));

end

W_power = mean(W_power,2);

if tf_data==1
    
    if tf_hemisphere == 1
    
    save('data5/asymmetry_erai_nh.mat','lambda_theory_av','lambda_era_av',...
       'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
    else
    
    save('data5/asymmetry_erai_sh.mat','lambda_theory_av','lambda_era_av',...
       'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
   
    end
    
else
    
      if tf_hemisphere ==1
    
    save('data5/asymmetry_era5_nh.mat','lambda_theory_av','lambda_era_av',...
       'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
    else
    
    save('data5/asymmetry_era5_sh.mat','lambda_theory_av','lambda_era_av',...
      'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
   
      end
    
end

% if tf_data==1
%     
%     if tf_hemisphere == 1
%     
%     save('data/asymmetry_erai_nh.mat','lambda_theory_av','lambda_era_av',...
%        'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
%     else
%     
%     save('data/asymmetry_erai_sh.mat','lambda_theory_av','lambda_era_av',...
%        'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
%    
%     end
%     
% else
%     
%       if tf_hemisphere ==1
%     
%     save('data/asymmetry_era5_nh.mat','lambda_theory_av','lambda_era_av',...
%        'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
%     else
%     
%     save('data/asymmetry_era5_sh.mat','lambda_theory_av','lambda_era_av',...
%        'lambda_era_mid_av','r_mid_av','L_mid_av','k_av','T_surf_av','W_power')
%    
%       end
%     
% end
        
    
% Make Plots of Seasonal Cycle
% figure(1)
% plot(lambda_era_av,'Linewidth',1.5); hold on;
% plot(lambda_theory_av,'linewidth',1.5)
% legend('ERA','Theory'); legend boxoff
% ylabel('\lambda')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% 
% 
% figure(2)
% plot(Wavenumber_av.*L_mid_av,'Linewidth',1.5)
% ylabel('k')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% 
% 
% figure(3)
% plot(r_mid_av,'Linewidth',1.5)
% ylabel('r')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% 
% 
% figure(4)
% plot(L_mid_av/1000,'Linewidth',1.5)
% ylabel('L_D(km)')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% 
% figure(5)
% plot(T_surf_av,'Linewidth',1.5)
% ylabel('T_{g}(K)')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
