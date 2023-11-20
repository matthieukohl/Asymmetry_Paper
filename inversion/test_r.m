%% Test File to check the lambda(r) behavior
clear;
close all;

list = {'rescale_true0.01','rescale0.02','rescale0.05','rescale0.1','rescale0.2','rescale0.3',...                     
'rescale0.4','rescale_true0.6','rescale0.8','rescale_true1.0'};

list1 = {'_day0225h00.nc','_day0250h00.nc','_day0275h00.nc','_day0300h00.nc'};

%list2 = {'_day0125h00.nc','_day0150h00.nc','_day0175h00.nc','_day0200h00.nc'};

%list3 = {'_day0325h00.nc','_day0350h00.nc','_day0375h00.nc','_day0400h00.nc'};




r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

Asymmetry = zeros(length(r),1);

% Constants

Omega = 7.2921e-5; % Rotation Rate
R = 6371000.0; % Planetary Radius 
g = 9.81; % gravitational constant
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

% Pre-allocate Stratospheric Level

strat_level = zeros(length(list),1);

%strat_level = [18;14;13;14;14;13;12;18;12;18];


for tt = 1:10
tt   

omega_full = [];

for i=1:4
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_',list{tt},'/temp',list1{i});
path_omega = strcat('/disk7/mkohl/interpolation/output1/data_',list{tt},'/omega',list1{i});

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
%event_latspan = [21:46];
event_lonspan = [1:256]; % 1:256
event_timespan = [1:100];
level_indices = [1:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);
Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

omega_in= ncread(path_omega,'omega_plevel', start, count);
temp_in= ncread(path_temp,'temp_plevel', start, count);


[omega,T] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        omega(:, :, k, t) = omega_in(:, :, k, t)';
        T(:, :, k, t) = temp_in(:, :, k, t)';
    end
end

%lapse_rate = moist_adiabatic_lapse_rate(T,Level,'simple');

T = mean(T(:,:,:,:),4);

% p_sat - saturation pressure water vapor

p_sat = p_triple*exp(-L/Rc*(1./T(:,:,:)-1/T_triple));

% r_sat - saturated mixing ratio

r_sat = epsilon*p_sat./(Level-p_sat);

% Calculate Virtual Temperature

T_v = T.*(1+0.608*r_sat);

% Calculate Lapse in height coordinates

%lapse_rate = -g*Level./(Ra*T_v).*d_dp(T,level);

lapse_rate = lapse_rate_function(T,Level,level);


% lapse_rate_mean = mean(mean(lapse_rate,1),2);
% indx = lapse_rate_mean<-2*10^(-3);
% strat_level_place(i) = max(find(indx));

% figure(2)
% plot(squeeze(mean(mean(lapse_rate(:,:,:),1),2)),level)
% set(gca, 'YDir','reverse')
% % 
% figure(3)
% plot(squeeze(mean(mean(T(:,:,:),1),2)),level)
% set(gca, 'YDir','reverse')

% Find the Stratosphere

for ii=1:length(lat)
    for jj=1:length(lon)
        indx = lapse_rate(ii,jj,:)<-2*10^(-3);
        max_level(ii,jj) = max(find(indx));
    end   
end

strat_level_place(i) = round(mean(mean(max_level)),0);

omega_full = cat(4,omega_full,omega(:,:,:,:));
a(i) = mean(Lambda_weighted(omega(:,:,1:strat_level_place(i),:),lat));
end 

strat_level(tt) = round(mean(strat_level_place),0);
        
% Calculate the Asymmetry Factor vs. r 

omega_full = omega_full - mean(omega_full(:,:,:,:),2);

Asymmetry(tt) = mean(Lambda(omega_full(:,:,10,:)));

end


figure(1)
semilogx(r,Asymmetry); hold on
title('Asymmetry vs. Abs for Omega-Fields')
xlabel('r')
ylabel('Lambda')
% savefig('mode_r_2.fig')


% Asymmetry = [0.7030,0.7030,0.7030,0.7030,0.6944,0.6656,0.6373,0.5861,0.5432,0.5122];
% 
% figure(2)
% plot(squeeze(mean(mean(lapse_rate(:,:,:),1),2)),level)
% set(gca, 'YDir','reverse')
% 
% figure(3)
% plot(squeeze(mean(mean(T(:,:,:),1),2)),level)
% set(gca, 'YDir','reverse')

strat_level = [18,18,18,18,17,17,16,16,16,16];


