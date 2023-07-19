%% Test File to check the lambda(r) behavior
clear; 
close all;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

abs = [0.4,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,4.0,6.0];

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808 ...
310.7297, 316.1383];

Asymmetry = zeros(length(list),1);

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


for tt = 1:15
tt    
for i=1:6

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

event_latspan = [93:108]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [1:200];
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

lapse_rate = -g*Level./(Ra*T_v).*d_dp(T,level);

% Find the Stratosphere

for ii=1:length(lat)
    for jj=1:length(lon)
        indx = lapse_rate(ii,jj,:)<-2*10^(-3);
        max_level(ii,jj) = max(find(indx));
    end 
end

strat_level_place(i) = round(mean(mean(max_level)),0);

end 

strat_level(tt) = round(mean(strat_level_place),0);
        
% Calculate the Asymmetry Factor vs. r 

omega = omega - mean(omega(:,:,:,:),2);

Asymmetry(tt) = mean(Lambda(omega(:,:,1:strat_level,:)));

end


% figure(1)
% plot(surf_temp,Asymmetry)
% title('Asymmetry vs. Abs for Omega-Fields')
% xlabel('abs')
% ylabel('Lambda')


%strat_level_mode =  [8; 9; 10; 10 ; 11; 11; 12; 12; 13; 13; 14; 15; 15; 16; 18];

%strat_level_turbulent = [8,9,10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15,16,18];

