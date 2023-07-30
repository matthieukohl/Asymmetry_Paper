% Load the power spectrum for omega from ncl calculation


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
L = 2493*10^3; % Latent Heat Water Vapour
g = 9.80665;  
epsilon = 0.621; % Mass ratio Water Vapour / Air

% file = ['/net/graupel/archive1/users/mkohl/reanalysis/asymmetry_project/erai/'...
% 'power_spectrum_omega_2015_ncep2.nc'];

% file = ['/net/graupel/archive1/users/mkohl/reanalysis/asymmetry_project/erai/'...
% 'power_spectrum_erai_omega_NH_2015.nc'];

file = ['/net/graupel/archive1/users/mkohl/reanalysis/asymmetry_project/erai/'...
'power_spectrum_era5_NH_6_2015.nc'];


lat = ncread(file,'latitude'); lat = double(lat); 
lon = ncread(file,'longitude');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; 

% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

time = ncread(file,'time');

[a,p200] = min(abs(level-200*1e2));
[a,p500] = min(abs(level-500*1e2));
[a,p1000] = min(abs(level-1000*1e2));

power_w = ncread(file,'power');
power_w = mean(power_w,3); % time average
power_w = mean(power_w(:,p500),2); % pressure average
power_w(1) = 0;


k = 1:1:length(lat);
k = k';

loglog(k,power_w); hold on;
xlabel('latitudinal wavenumber')
ylabel('power')

% convert to total zonal wavenumber

k_centroid = sum(k.*power_w)/sum(power_w);

n = k_centroid-1;

K = sqrt((n+1)*n/R^2);
%Ld = 1e6/(2*sqrt(2));
Ld = 5.8*1e5;

K = K*Ld