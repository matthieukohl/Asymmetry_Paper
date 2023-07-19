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

lat = ncread(path_temp,'latitude');
lon = ncread(path_temp,'longitude');
time = ncread(path_temp,'time');

years = 2009:1:2009;

season = 'jja';  

counter = 0;

k_rhs = 0;
skew_summer = 0;
spec_rhs = zeros(1,length(lon));

for tt = 1:length(years)
  

timespan = season_selector(path_temp, season, years(tt), years(tt));


dphi = (lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = (lon(2) - lon(1)) / 180.0 * 3.1415926;

lat = flip(lat);

[Lon,Lat] = meshgrid(lon,lat);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [1, 1, timespan(1)];
count = [length(lon), length(lat),length(timespan)];

T_in= ncread(path_temp,'t', start, count);
u_in= ncread(path_u,'u', start, count);
v_in= ncread(path_v,'v', start, count);


T = rearrange(T_in,lat,lon,1,timespan);
u = rearrange(u_in,lat,lon,1,timespan);
v = rearrange(v_in,lat,lon,1,timespan);

  
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


skew_summer_new = Skewness(rhs(:,:,t));

%[k_rhs_new,spec_rhs_new] = Spectrum1(rhs(:,:,1),lat,lon);
[k_rhs_new,spec_rhs_new] = Spectrum2(rhs(:,:,t),lat,lon);

k_rhs = k_rhs + k_rhs_new;
spec_rhs = spec_rhs + spec_rhs_new;
skew_summer = skew_summer + skew_summer_new;
counter = counter + 1;
end

end


k_rhs = k_rhs/counter;
spec_rhs = spec_rhs/counter;
k_rhs = 0.5*k_rhs*1e6/sqrt(2);

k_rhs_summer = k_rhs;
spec_rhs_summer = spec_rhs;

skew_summer = skew_summer/counter;

% winter

counter = 0;

k_rhs = 0;
skew_winter = 0;
spec_rhs = zeros(1,length(lon));

%list_month = {'jan','feb','dec'};
list_month = {'ann'};


for tt = 1:length(years)
    
for ii = 1:length(list_month)
    
season = list_month{ii}; 
  

timespan = season_selector(path_temp, season, years(tt), years(tt));


dphi = (lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = (lon(2) - lon(1)) / 180.0 * 3.1415926;

lat = flip(lat);

[Lon,Lat] = meshgrid(lon,lat);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [1, 1, timespan(1)];
count = [length(lon), length(lat),length(timespan)];

T_in= ncread(path_temp,'t', start, count);
u_in= ncread(path_u,'u', start, count);
v_in= ncread(path_v,'v', start, count);

T = rearrange(T_in,lat,lon,1,timespan);
u = rearrange(u_in,lat,lon,1,timespan);
v = rearrange(v_in,lat,lon,1,timespan);

  
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
[k_rhs_new,spec_rhs_new] = Spectrum2(rhs(:,:,t),lat,lon);

k_rhs = k_rhs + k_rhs_new;
spec_rhs = spec_rhs + spec_rhs_new;
skew_winter = skew_winter + skew_winter_new;
counter = counter + 1;
end

end


end



R_eff = cosd(45)*R;
N = length(lon);
%k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];

k = 0.5*k*1e6/sqrt(2);

k_rhs = k_rhs/counter;
spec_rhs = spec_rhs/counter;
k_rhs = 0.5*k_rhs*1e6/sqrt(2);

k_rhs_winter = k_rhs;
spec_rhs_winter = spec_rhs;

skew_winter = skew_winter/counter;


figure(1)
x0=10; y0=10; width=1000; height=600; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)
loglog(k,spec_rhs_winter,'b'); hold on;
loglog(k,spec_rhs_summer,'r'); hold on;
loglog(k,1e-25*k.^(-3/2),'k');
legend('winter','summer','k^{-3/2}'); legend boxoff
xlabel('k')

indL = 800;
indR = 1300;

subplot(1,2,2)
loglog(k(indL:indR),spec_rhs_winter(indL:indR),'b'); hold on;
loglog(k(indL:indR),spec_rhs_summer(indL:indR),'r'); hold on;
loglog(k(indL:indR),1e-25*k(indL:indR).^(-3/2),'k');
legend('winter','summer','k^{-3/2}'); legend boxoff
xlabel('k')

figure(2)
semilogx(k,spec_rhs_winter/max(spec_rhs_winter)); hold on
%semilogx(k,spec_rhs_summer/max(spec_rhs_summer));
%legend('winter','summer'); legend boxoff
xlabel('k')
title('\rm spec rhs era5')

indL = 30;
indR = 170;

figure(3)
loglog(k,spec_rhs_winter/max(spec_rhs_winter)); hold on
loglog(k(indL:indR),0.5*k(indL:indR).^1/4);
%semilogx(k,spec_rhs_summer/max(spec_rhs_summer));
%legend('winter','summer'); legend boxoff
xlabel('k')
title('\rm spec rhs era5')
legend('era5','k^{1/4}')
legend boxoff

save('spec_rhs_era5.mat','k','spec_rhs_winter');


saveas(gcf,'figures_test/spec_rhs_reanalysis','epsc')







