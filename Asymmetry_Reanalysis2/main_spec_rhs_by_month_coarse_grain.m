% Calculate the spectrum of the RHS for ERA5

clear; close all;

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
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% Define Path 
path_temp = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/temp_2009_2018_NH.nc';
path_u = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/u_2009_2018_NH.nc';
path_v = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/v_2009_2018_NH.nc';

% Read in Grid-Variables Lat/Lon 

lat_original = ncread(path_temp,'latitude');
lon_original = ncread(path_temp,'longitude');
time = ncread(path_temp,'time');
lat_original = flip(lat_original);

years = 2009:1:2018;

spec_rhs_era5 = zeros(240,12);

list_month = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};

for ii = 1:length(list_month)
    
counter = 0;

k_rhs = 0;
spec_rhs = zeros(1,240);
skew_winter = 0;  
    
for tt = 1:length(years)
    

    
season = list_month{ii}; 
  

timespan = season_selector(path_temp, season, years(tt), years(tt));

% Read in Variables (Lat/Lon need to be flipped)

start = [1, 1, timespan(1)];
count = [length(lon_original), length(lat_original),length(timespan)];

T_in= ncread(path_temp,'t', start, count);
u_in= ncread(path_u,'u', start, count);
v_in= ncread(path_v,'v', start, count);

T = rearrange(T_in,lat_original,lon_original,1,timespan);
u = rearrange(u_in,lat_original,lon_original,1,timespan);
v = rearrange(v_in,lat_original,lon_original,1,timespan);


% coarse grain

u = movmean(movmean(u,7,1),7,2); % averaging

u = u(1:6:end,1:6:end,:); % downsampling

v = movmean(movmean(v,7,1),7,2); % averaging

v = v(1:6:end,1:6:end,:); % downsampling

T = movmean(movmean(T,7,1),7,2); % averaging

T = T(1:6:end,1:6:end,:); % downsampling

lat = lat_original(1:6:end);
lon = lon_original(1:6:end);

dphi = abs((lat(2) - lat(1))) / 180.0 * 3.1415926;
dlambda = abs((lon(2) - lon(1))) / 180.0 * 3.1415926;

[Lon,Lat] = meshgrid(lon,lat);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);



  
f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);
            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(lat), length(lon), length(timespan)));
    

[rhs] = deal(zeros(length(lat), length(lon),length(timespan)));


for t = 1 : length(timespan)
    
% term A (Q vector)


    
    du_dlambda(:, :, t) = d_dlambda(u(:, :, t), dlambda);
    dv_dlambda(:, :, t) = d_dlambda(v(:, :, t), dlambda);
    du_dphi(:, :, t) = d_dphi(u(:, :,t), dphi);
    dv_dphi(:, :, t) = d_dphi(v(:, :, t), dphi);

    dT_dlambda(:, :, t) = d_dlambda(T(:, :, t), dlambda);
    dT_dphi(:, :, t) = d_dphi(T(:, :, t), dphi);
    Qx(:, :, t) = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda(:, :, t) .* dT_dlambda(:, :,t) + ...
        (1 ./ cos(Phi))    .* dv_dlambda(:, :, t) .* dT_dphi(:, :, t));
    Qy(:, :, t) = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi(:, :, t)    .* dT_dlambda(:, :, t) + ...
                              dv_dphi(:, :, t)    .* dT_dphi(:, :, t));
    div_temp = div(Qx(:, :,t), Qy(:, :,t), Phi, dphi, dlambda);
    A(:, :, t) = 2 * Ra / (500*1e2) * div_temp;
    
    
clear('temp_x', 'temp_y', 'div_temp');


% right hand side

rhs(:, :, t) = A(:, :, t); 

skew_winter_new = Skewness(rhs(:,:,t));


%[k_rhs_new,spec_rhs_new] = Spectrum1(rhs(:,:,1),lat,lon);
[k_rhs_new,spec_rhs_new] = Spectrum1_weighted(rhs(:,:,t),lat,lon);
%[k_rhs_new,spec_rhs_new] = Spectrum2(rhs(:,:,t),lat,lon);

k_rhs = k_rhs + k_rhs_new;
spec_rhs = spec_rhs + spec_rhs_new;
skew_winter = skew_winter + skew_winter_new;
counter = counter + 1;
end

end

spec_rhs_era5(:,ii) = spec_rhs'/counter;
spec_rhs_era5(:,ii) = spec_rhs_era5(:,ii)/max(spec_rhs_era5(:,ii),[],1);


end

save('spec_rhs_NH_era5_coarse_1p5_res.mat','lat','lon','spec_rhs_era5');

