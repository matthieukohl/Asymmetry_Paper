% Calculate the reduction factor and predicted asymmetry
% Based on x average for era5

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

% Input for Era5 files

years  = [2009:2009];
%years = [2014:2014];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'north';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[wavenumber_total_vs_month,wavenumber_up_vs_month,wavenumber_down_vs_month] = deal(zeros(12,1));
else
[wavenumber_total_vs_month,wavenumber_up_vs_month,wavenumber_down_vs_month] = deal(zeros(12,1));
end

tic
for ii = 1:length(months)
ii
season = months{ii};

counter = 0;

wavenumber_total = 0;
wavenumber_up = 0;
wavenumber_down = 0;

for year = years

disp('year:')
disp(year)


% file stored in aimsir
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/ncep2/omega.', num2str(year),'.nc'];


lat = ncread(file,'lat'); lat = double(lat); 
lon = ncread(file,'lon');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; 

[a,lat70] = min(abs(lat-70));
[a,lat30] = min(abs(lat-30));

lat = lat(lat70:lat30);

lat = flip(lat);
level = flip(level);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);
dphi = abs(lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = abs(lon(2) - lon(1)) / 180.0 * 3.1415926;

% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

tic

for t = time_range

start = [1, lat70, 1, t];
count = [length(lon), abs(lat70-lat30)+1, length(level), 1];

% have format lon/lat/level/time

omega_in = ncread(file,'omega',start,count);

omega = rearrange_era(omega_in,lat,lon,level);

clear('omega_in','T_in');

% calculate the laplacian

[a,p500] = min(abs(level-500*1e2));

% now evaluate mean(lap(w)^2)/mean(w^2)
% over the whole domain/ over ascending regions
% over descending regions


w_2d = omega(:,:,p500);        
lap_w = spherical_laplacian(w_2d, phi, dphi, dlambda);


% paul's formulation

% square formulation

wavenumber_total = wavenumber_total + mean(mean(lap_w.^2))/(mean(mean(w_2d.^2)));
wavenumber_up = wavenumber_up + mean(mean(lap_w(w_2d<0).^2))/(mean(mean(w_2d(w_2d<0).^2)));
wavenumber_down = wavenumber_down + mean(mean(lap_w(w_2d>0).^2))/(mean(mean(w_2d(w_2d>0).^2)));

% absolute value formulation

% wavenumber_total = wavenumber_total + mean(mean(abs(lap_w)))/(mean(mean(abs(w_2d))));
% wavenumber_up = wavenumber_up + mean(mean(abs(lap_w(w_2d<0))))/(mean(mean(abs(w_2d(w_2d<0)))));
% wavenumber_down = wavenumber_down + mean(mean(abs(lap_w(w_2d>0))))/(mean(mean(abs(w_2d(w_2d>0)))));

% wavenumber_total = wavenumber_total + mean(mean(abs(w_2d)))/(mean(mean(abs(lap_w))));
% wavenumber_up = wavenumber_up + mean(mean(abs(w_2d(w_2d<0))))/(mean(mean(abs(lap_w(w_2d<0)))));
% wavenumber_down = wavenumber_down + mean(mean(abs(w_2d(w_2d>0))))/(mean(mean(abs(lap_w(w_2d>0)))));


% my formulation

% % absolute value
% 
%ratio = lap_w./w_2d;
% 
% wavenumber_total = wavenumber_total + mean(mean(abs(ratio)));
% wavenumber_up = wavenumber_up + mean(mean(abs(ratio(w_2d<0))));
% wavenumber_down = wavenumber_down + mean(mean(abs(ratio(w_2d>0))));

counter = counter + 1;

end
end

wavenumber_total = wavenumber_total/counter;
wavenumber_up = wavenumber_up/counter;
wavenumber_down = wavenumber_down/counter;

wavenumber_total_vs_month(ii) = (wavenumber_total/4)^(1/4);
wavenumber_up_vs_month(ii) = (wavenumber_up/4)^(1/4);
wavenumber_down_vs_month(ii) = (wavenumber_down/4)^(1/4);


end % end loop over time

wavenumber_total_estimate_vs_month = 2*wavenumber_up_vs_month.*wavenumber_down_vs_month./(wavenumber_up_vs_month+wavenumber_down_vs_month);

% nondimensionalize

load('data_theory/theory_ncep2_NH_april29th.mat'); 
load('data_theory/theory_ncep2_SH_april29th.mat');

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

wavenumber_total_vs_month = wavenumber_total_vs_month.*LD_NH_500;
wavenumber_total_estimate_vs_month = wavenumber_total_estimate_vs_month.*LD_NH_500;
wavenumber_down_vs_month = wavenumber_down_vs_month.*LD_NH_500; 
wavenumber_up_vs_month = wavenumber_up_vs_month.*LD_NH_500; 


figure
plot(wavenumber_down_vs_month,'b','Linewidth',1.4); hold on;
plot(wavenumber_up_vs_month,'r','Linewidth',1.4); hold on;
plot(wavenumber_total_vs_month,'k','Linewidth',1.4); hold on;
plot(wavenumber_total_estimate_vs_month,'k--','Linewidth',1.4);
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('Reduction factor r')
ylabel('Wavenumber k')
title('\rm mean((Lap w)^2) / mean(w^2)')
set(gca,'FontSize',12);
legend('k_d','k_u','k','k estimate'); legend boxoff;