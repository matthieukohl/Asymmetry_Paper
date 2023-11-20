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
[w_up_NH,S_NH,Ld_NH] = deal(zeros(37,12));
else
[w_up_SH,S_SH,Ld_SH] = deal(zeros(37,12));
end

lengths_vs_month = zeros(12,1);

tic
for ii = 1:1
ii
season = months{ii};

norm = 0;
counter_w_up = 0;
w_up = zeros(37,1);
S = zeros(37,1);

lengths = [];

for year = years

disp('year:')
disp(year)

if tf_hemisphere ==1
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_',num2str(ii),'_', num2str(year),'.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ii),'_',num2str(year),'.nc'];
end

lat = ncread(file,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(file,'longitude');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; level = flip(level);

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
T_in = ncread(file,'t',start,count);

T = rearrange_era(T_in,lat,lon,level);
omega = rearrange_era(omega_in,lat,lon,level);

clear('omega_in','T_in');

Level = repmat(reshape(level, 1, 1, length(level)), length(lat), length(lon), 1);
theta = T .* (p0 ./ Level).^ kappa;
sigma_accu = -Ra * T./ (Level .* theta) .* d_dp(theta, level);

r = r_parameter1(T,Level,level);
density = rho(T,Level);
lapse = -g*density.*d_dp(T,level)*10^3; % K/km

% find various levels

[~,p500] = min(abs(level-50000));
[~,p400] = min(abs(level-40000));
[~,p700] = min(abs(level-70000));
[~,plow] = min(abs(level-92500));
surf = 1;



norm = norm+1;

w = - omega(:,:,p500);

for tt = 1: size(w,1)

w_slice = w(tt,:);
w_slice(w_slice>0) = 0;
w_down = w_slice<0;

input_array = w_down; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths = cat(2,lengths,af-a0);
clear('input_array');

end

lengths_up = [];

for tt = 1: size(w,1)

w_slice = w(tt,:);
w_slice(w_slice<0) = 0;
w_up = w_slice>0;

input_array = w_up; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths_up = cat(2,lengths_up,af-a0);
clear('input_array');

end


end % end loop over time

toc

end % end loop over years

lengths_vs_month(ii) = mean(lengths);
length_vs_month_alternative(ii) = mean(1./lengths);
lengths_vs_month_up(ii) = mean(lengths_up);
lengths_vs_month_up_alternative(ii) = mean(1./lengths_up);

end

% convert index length to km
dphi = abs(lon(2)-lon(1));
lengths = lengths*2*pi*R*cosd(45)*dphi/360;
lengths_vs_month = lengths_vs_month*2*pi*R*cosd(45)*dphi/360;
lengths_vs_month_up = lengths_vs_month_up*2*pi*R*cosd(45)*dphi/360;

lengths_vs_month_alternatve = lengths_vs_month_alternative/(2*pi*R*cosd(45)*dphi/360);
lengths_vs_month_up_alternatve = lengths_vs_month_up_alternative/(2*pi*R*cosd(45)*dphi/360);


figure
histogram(lengths)

load('data_theory/theory_NH_dec7th.mat'); 
load('data_theory/theory_SH_dec7th.mat');

dp = 800*1e2;
f = 1e-4;

% LD_NH_500  = 1e6/(2*sqrt(2))*ones(12,1);
% LD_SH_500 = LD_NH_500;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % 

Ld = mean(LD_NH_500);
lengths = lengths/Ld;
lengths_vs_month = lengths_vs_month./LD_NH_500;
lengths_vs_month_up = lengths_vs_month_up./LD_NH_500;


figure
histogram(lengths,'normalization','pdf')
set(gca,'YScale','log')
xlabel('Ld')
title(['Ld=',num2str(round(mean(lengths),2)),' kd=',num2str(round(1/mean(lengths),2)),'(',num2str(round(mean(1./lengths),2)),') ERA5',])
set(gca,'FontSize',12)



lengths = 1./lengths;


lengths_vs_month = 1./lengths_vs_month;



figure
plot(lengths_vs_month)
