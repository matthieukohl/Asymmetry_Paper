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

years  = [2009:2011];
%years = [2014:2014];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'NH';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[wavenumber_total_vs_month,wavenumber_up_vs_month,wavenumber_down_vs_month,length_ascent_vs_month,length_descent_vs_month] = deal(zeros(12,1));
else
[wavenumber_total_vs_month,wavenumber_up_vs_month,wavenumber_down_vs_month,length_ascent_vs_month,length_descent_vs_month] = deal(zeros(12,1));
end

tic
for ii = 1:length(months)
ii
season = months{ii};

counter = 0;

wavenumber_total = 0;
wavenumber_up = 0;
wavenumber_down = 0;

lengths_down = [];
lengths_up = [];

for year = years

disp('year:')
disp(year)

% % era5
% if tf_hemisphere ==1
% file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_',num2str(ii),'_', num2str(year),'.nc'];
% else
% file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ii),'_',num2str(year),'.nc'];
% end

% erai

if tf_hemisphere ==1
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai_test/erai_NH_omega_', num2str(year),'.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/erai_test/erai_SH_omega_', num2str(year),'.nc'];
end

lat = ncread(file,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(file,'longitude');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; level = flip(level);


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

start = [1, 1, 1, t];
count = [length(lon), length(lat), length(level), 1];

start_sf = [1, 1, t];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega_in = ncread(file,'w',start,count);

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

% start calculating

% [a,p500] = min(abs(level-500e2));
%     
% w = - squeeze(omega(:,:,p500));

% for i = 1: size(w,1)
% 
% w_slice = squeeze(w(i,:));
% %w_slice(w_slice>0) = 0;
% w_down = w_slice<0;
% 
% input_array = w_down; 
% input_array = [0 input_array 0];
% 
% input_diff = diff(input_array);
% 
% a0 = find(input_diff==1);
% af = find(input_diff==-1);
% 
% 
% %lengths_down = cat(2,lengths_down,(af-a0));
% 
% index_length =  2*pi*R*cosd(lat(i))*(lon(2)-lon(1))/360;
% 
% lengths_down = cat(2,lengths_down,(af-a0)*index_length);
% clear('input_array','a0','af');
% 
% end
% 
% for i = 1: size(w,1)
% 
% w_slice = squeeze(w(i,:));
% %w_slice(w_slice>0) = 0;
% w_up = w_slice>0;
% 
% input_array = w_up; 
% input_array = [0 input_array 0];
% 
% input_diff = diff(input_array);
% 
% a0 = find(input_diff==1);
% af = find(input_diff==-1);
% 
% %lengths_up = cat(2,lengths_up,(af-a0));
% 
% index_length =  2*pi*R*cosd(lat(i))*(lon(2)-lon(1))/360;
% 
% lengths_up = cat(2,lengths_up,(af-a0)*index_length);
% 
% clear('input_array','a0','af');
% 
% end


counter = counter + 1;

end


end


wavenumber_total = wavenumber_total/counter;
wavenumber_up = wavenumber_up/counter;
wavenumber_down = wavenumber_down/counter;

wavenumber_total_vs_month(ii) = (wavenumber_total/4)^(1/4);
wavenumber_up_vs_month(ii) = (wavenumber_up/4)^(1/4);
wavenumber_down_vs_month(ii) = (wavenumber_down/4)^(1/4);

% wavenumber_total_vs_month(ii) = (wavenumber_total/2)^(1/2);
% wavenumber_up_vs_month(ii) = (wavenumber_up/2)^(1/2);
% wavenumber_down_vs_month(ii) = (wavenumber_down/2)^(1/2);



% length_descent_vs_month(ii) = mean(lengths_down);
% length_ascent_vs_month(ii) = mean(lengths_up);

end % end loop over time

wavenumber_total_estimate_vs_month = 2*wavenumber_up_vs_month.*wavenumber_down_vs_month./(wavenumber_up_vs_month+wavenumber_down_vs_month);

% nondimensionalize

load('data_theory/theory_NH_dec7th.mat'); 
load('data_theory/theory_SH_dec7th.mat');

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

% index_length =  2*pi*R*cosd(45)*(lon(2)-lon(1))/360;
% 
% %length_descent_vs_month = length_descent_vs_month*index_length;
% length_descent_vs_month = length_descent_vs_month./LD_NH_500;
% 
% %length_ascent_vs_month = length_ascent_vs_month*index_length;
% length_ascent_vs_month = length_ascent_vs_month./LD_NH_500;


figure
plot(wavenumber_down_vs_month,'b','Linewidth',1.4); hold on;
plot(wavenumber_up_vs_month,'r','Linewidth',1.4); hold on;
plot(wavenumber_total_vs_month,'k','Linewidth',1.4); hold on;
plot(wavenumber_total_estimate_vs_month,'k--','Linewidth',1.4);
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('Reduction factor r')
ylabel('Wavenumber k')
%title('\rm mean((Lap w)^2) / mean(w^2)')
set(gca,'FontSize',12);
legend('k_d','k_u','k','k estimate'); legend boxoff;

figure

plot(1./length_descent_vs_month,'b','Linewidth',1.4); hold on;
plot(1./length_ascent_vs_month,'r','Linewidth',1.4); hold on;
plot(2*(1./length_ascent_vs_month).*(1./length_descent_vs_month)./(1./length_ascent_vs_month+1./length_descent_vs_month),'k');
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('Reduction factor r')
ylabel('Wavenumber k')
set(gca,'FontSize',12);
legend('k_d','k_u','k'); legend boxoff;
title('\rm ERAI')
