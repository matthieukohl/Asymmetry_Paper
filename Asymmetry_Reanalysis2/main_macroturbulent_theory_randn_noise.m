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

years      = [2014:2018];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'south';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[lambda_NH_theory_mid,lambda_NH_theory_500,lambda_NH_theory_av,...
k_NH_mid,k_NH_500,k_NH_av,L_mid_NH,L_500_NH,L_av_SH,...
r_mid_NH,r_500_NH,r_av_NH,trop_NH] = deal(zeros(12,1));
[spec_NH_mid,spec_NH_500,spec_NH_av] = deal(zeros(1440,12));
else
[lambda_SH_theory_mid,lambda_SH_theory_500,lambda_SH_theory_av,...
k_SH_mid,k_SH_500,k_SH_av,L_mid_SH,L_500_SH,L_av_SH,...
r_mid_SH,r_500_SH,r_av_SH,trop_SH] = deal(zeros(12,1));
[spec_SH_mid,spec_SH_500,spec_SH_av] = deal(zeros(1440,12));
end

tic
for ii = 1:length(months)

season = months{ii};

norm = 0;
counter = 0;
lambda_mid = 0; lambda_av = 0; lambda_500 = 0;
k_mid = 0; k_av = 0; k_500 = 0;
spec_mid = 0; spec_av = 0; spec_500 = 0;
r_mid = 0; r_av = 0; r_500 = 0;
L_mid = 0; L_av = 0; L_500 = 0;
trop = 0;


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
t
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

% find troposphere based on 2 K/km criteria

lapse_mean = mean(mean(lapse,1),2);

trop_new = find(lapse_mean>-2,1);

trop = trop+trop_new;

clear('lapse','theta','density','Level') 

% find various levels

[~,p500] = min(abs(level-50000));
[~,plow] = min(abs(level-92500));
surf = 1;
[~,pmid] = min(abs(level-(level(surf)-(level(surf)-level(trop_new))/2)));

% find the dominant Wavenumber

[k_av_new,spec_av_new] = Spectrum1(omega(:,:,plow:trop_new),lat,lon);
[k_mid_new,spec_mid_new] = Spectrum1(omega(:,:,pmid),lat,lon);
[k_500_new,spec_500_new] = Spectrum1(omega(:,:,p500),lat,lon);

k_av = k_av+k_av_new; spec_av = spec_av + spec_av_new;
k_mid = k_mid+k_mid_new; spec_mid = spec_mid + spec_mid_new;
k_500 = k_500+k_500_new; spec_500 = spec_500 + spec_500_new;

% find the Non-Dimensionalization to use

[L_mid_new,L_500_new,L_av_new] = Deformation_Radius(sigma_accu,f,level,lat,trop_new);

L_mid = L_mid+L_mid_new;

L_500 = L_500+L_500_new;

L_av = L_av + L_av_new;

clear('sigma_accu')

% find the r-value to use

[r_mid_new,r_500_new,r_av_new] = reduction_factor(r,level,trop_new);

r_mid = r_mid+r_mid_new;

r_500 = r_500+r_500_new;

r_av = r_av + r_av_new;

clear('r')

% calculate theoretical prediction

% do multiple iterations

for jj = 1:20

lambda_mid_new = toy_model_randn(r_mid_new,0.5);

lambda_500_new = toy_model_randn(r_500_new,0.5);

lambda_av_new = toy_model_randn(r_av_new,0.5);

lambda_mid = lambda_mid+lambda_mid_new;

lambda_500 = lambda_500 + lambda_500_new;

lambda_av = lambda_av + lambda_av_new;

counter = counter+1;

end

norm = norm + 1;

%lambda_test = lambda_test+Lambda_x(omega);

end % end loop over time

toc

end % end loop over years



lambda_mid = squeeze(lambda_mid/counter); lambda_500 = squeeze(lambda_500/counter);
lambda_av = squeeze(lambda_av/counter);
k_mid = k_mid/norm; k_500 = k_500/norm; k_av = k_av/norm;
spec_mid = squeeze(spec_mid)/norm; spec_500 = squeeze(spec_500)/norm;
spec_av = squeeze(spec_av)/norm;
L_mid = L_mid/norm; L_500 = L_500/norm; L_av = L_av/norm;
r_mid = r_mid/norm; r_500 = r_500/norm; r_av = r_av/norm;
trop = trop/norm;


if tf_hemisphere ==1
lambda_NH_theory_mid(ii) = lambda_mid;
lambda_NH_theory_500(ii) = lambda_500;
lambda_NH_theory_av(ii) = lambda_av;
k_NH_mid(ii) = k_mid; k_NH_500(ii) = k_500; k_NH_av(ii) = k_av;
spec_NH_mid(:,ii) = spec_mid; spec_NH_500(:,ii) = spec_500;
spec_NH_av(:,ii) = spec_av;
L_mid_NH(ii) = L_mid; L_500_NH(ii) = L_500; L_av_NH(ii) = L_av;
r_mid_NH(ii) = r_mid; r_500_NH(ii) = r_500; r_av_NH(ii) = r_av;
trop_NH(ii) = trop;
else
lambda_SH_theory_mid(ii) = lambda_mid;
lambda_SH_theory_500(ii) = lambda_500;
lambda_SH_theory_av(ii) = lambda_av;
k_SH_mid(ii) = k_mid; k_SH_500(ii) = k_500; k_SH_av(ii) = k_av;
spec_SH_mid(:,ii) = spec_mid; spec_SH_500(:,ii) = spec_500;
spec_SH_av(:,ii) = spec_av;
L_mid_SH(ii) = L_mid; L_500_SH(ii) = L_500; L_av_SH(ii) = L_av;
r_mid_SH(ii) = r_mid; r_500_SH(ii) = r_500; r_av_SH(ii) = r_av;
trop_SH(ii) = trop;
end


end
toc

if tf_hemisphere ==1

 save('data_theory/theory_NH_rand_alpha0.5.mat', ...
      'lambda_NH_theory_mid','lambda_NH_theory_500','lambda_NH_theory_av',...
      'k_NH_mid','k_NH_av','k_NH_500','spec_NH_mid','spec_NH_av','spec_NH_500',...
  'L_mid_NH','L_av_NH','L_500_NH','r_mid_NH','r_500_NH','r_av_NH',...
  'trop_NH', 'lat', 'level')
  
else
   save('data_theory/theory_SH_rand_alpha0.5.mat', ...
      'lambda_SH_theory_mid','lambda_SH_theory_500','lambda_SH_theory_av',...
      'k_SH_mid','k_SH_av','k_SH_500','spec_SH_mid','spec_SH_av','spec_SH_500',...
  'L_mid_SH','L_av_SH','L_500_SH','r_mid_SH','r_500_SH','r_av_SH',...
  'trop_SH', 'lat', 'level')
    
end