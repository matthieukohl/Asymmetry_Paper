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

years      = [2009:2018];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'north';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[lambda_NH_theory_mid,lambda_NH_theory_500,lambda_NH_theory_av,...
k_NH_mid,k_NH_500,k_NH_av,L_mid_NH,L_500_NH,L_av_NH,...
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
file_u = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/u_2009_2018_NH.nc'];
file_v = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/v_2009_2018_NH.nc'];
file_temp = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/temp_2009_2018_NH.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ii),'_',num2str(year),'.nc'];
file_u = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/u_2009_2018_SH.nc'];
file_v = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/v_2009_2018_SH.nc'];
file_temp = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/temp_2009_2018_SH.nc'];
end

lat = ncread(file_u,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(file_u,'longitude');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; level = flip(level);

dphi = abs((lat(2) - lat(1))) / 180.0 * 3.1415926;
dlambda = abs((lon(2) - lon(1))) / 180.0 * 3.1415926;

[Lon,Lat] = meshgrid(lon,lat);

phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);


% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

time = ncread(file,'time');

time_range = season_selector(file_u, season, year, year);


for t = 1:length(time_range)

start = [1, 1, time_range(t)];
count = [length(lon), length(lat),1];

start_sf = [1, 1, 1, t];
count_sf = [length(lon), length(lat), length(level),1];

% have format lon/lat/level/time
u500_in = ncread(file_u,'u',start,count);
v500_in = ncread(file_v,'v',start,count);
T500_in = ncread(file_temp,'t',start,count);
T_in = ncread(file,'t',start_sf,count_sf);

u500 = rearrange(u500_in,lat,lon,level,1);
v500 = rearrange(v500_in,lat,lon,level,1);
T500 = rearrange(T500_in,lat,lon,level,1);
T = rearrange_era(T_in,lat,lon,level);

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

clear('T','lapse','theta','density','Level') 


% calculate the rhs

% term A (Q vector)

    du_dlambda = d_dlambda(u500, dlambda);
    dv_dlambda = d_dlambda(v500, dlambda);
    du_dphi = d_dphi(u500, dphi);
    dv_dphi = d_dphi(v500, dphi);

    dT_dlambda = d_dlambda(T500, dlambda);
    dT_dphi = d_dphi(T500, dphi);
    Qx = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda.* dT_dlambda + ...
        (1 ./ cos(Phi))    .* dv_dlambda.* dT_dphi);
    Qy = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi    .* dT_dlambda + ...
                              dv_dphi   .* dT_dphi);
    div_temp = div(Qx, Qy, Phi, dphi, dlambda);
    A = 2 * Ra / (500*1e2) * div_temp;
    
    
clear('temp_x', 'temp_y', 'div_temp');


% right hand side

rhs = A; 

% find various levels

[~,p500] = min(abs(level-50000));
[~,plow] = min(abs(level-92500));
surf = 1;
[~,pmid] = min(abs(level-(level(surf)-(level(surf)-level(trop_new))/2)));

% find the dominant Wavenumber


[k_500_new,spec_500_new] = Spectrum1(rhs,lat,lon);

k_500 = k_500+k_500_new; spec_500 = spec_500 + spec_500_new;

% find the Non-Dimensionalization to use

[L_mid_new,L_500_new,L_av_new] = Deformation_Radius(sigma_accu,f,level,lat,trop_new);

L_mid = L_mid + L_mid_new;
L_500 = L_500+L_500_new;
L_av = L_av + L_av_new;


clear('sigma_accu')

% find the r-value to use

[r_mid_new,r_500_new,r_av_new] = reduction_factor(r,level,trop_new);

r_500 = r_500+r_500_new;
r_mid = r_mid + r_mid_new;
r_av = r_av + r_av_new;


clear('r')

% calculate theoretical prediction


lambda_500_new = toy_model(r_500_new,k_500_new,L_500_new);
lambda_mid_new = toy_model(r_mid_new,k_500_new,L_500_new);
lambda_av_new = toy_model(r_av_new,k_500_new,L_500_new);


lambda_500 = lambda_500 + lambda_500_new;
lambda_mid = lambda_mid + lambda_mid_new;
lambda_av = lambda_av + lambda_av_new;

norm = norm+1;

%lambda_test = lambda_test+Lambda_x(omega);

end % end loop over time

toc

end % end loop over years



lambda_500 = squeeze(lambda_500/norm);
lambda_mid = squeeze(lambda_mid/norm);
lambda_av = squeeze(lambda_av/norm);

k_500 = k_500/norm; 
spec_500 = squeeze(spec_500)/norm;

L_500 = L_500/norm; 
L_mid = L_mid/norm;
L_av = L_av/norm;
r_500 = r_500/norm; 
r_mid = r_mid/norm;
r_av = r_av/norm;
trop = trop/norm;


if tf_hemisphere ==1
lambda_NH_theory_500(ii) = lambda_500;
lambda_NH_theory_mid(ii) = lambda_mid;
lambda_NH_theory_av(ii) = lambda_av;

k_NH_500(ii) = k_500; 
spec_NH_500(:,ii) = spec_500;

L_500_NH(ii) = L_500; 
L_mid_NH(ii) = L_mid;
L_av_NH(ii) = L_av;
r_500_NH(ii) = r_500; 
r_mid_NH(ii) = r_mid;
r_av_NH(ii) = r_av;
trop_NH(ii) = trop;
else

lambda_SH_theory_500(ii) = lambda_500;
lambda_SH_theory_mid(ii) = lambda_mid;
lambda_SH_theory_av(ii) = lambda_av;

k_SH_500(ii) = k_500; 
spec_SH_500(:,ii) = spec_500;

L_500_SH(ii) = L_500;
L_mid_SH(ii) = L_mid;
L_av_SH(ii) = L_av;

r_500_SH(ii) = r_500;
r_mid_SH(ii) = r_mid;
r_av_SH(ii) = r_av;
trop_SH(ii) = trop;
end


end
toc

if tf_hemisphere ==1

 save('data_theory/theory_NH_sin_k_rhs.mat', ...
      'lambda_NH_theory_mid','lambda_NH_theory_500','lambda_NH_theory_av',...
      'k_NH_mid','k_NH_av','k_NH_500','spec_NH_mid','spec_NH_av','spec_NH_500',...
  'L_mid_NH','L_av_NH','L_500_NH','r_mid_NH','r_500_NH','r_av_NH',...
  'trop_NH', 'lat', 'level')
  
else
   save('data_theory/theory_SH_sin_k_rhs.mat', ...
      'lambda_SH_theory_mid','lambda_SH_theory_500','lambda_SH_theory_av',...
      'k_SH_mid','k_SH_av','k_SH_500','spec_SH_mid','spec_SH_av','spec_SH_500',...
  'L_mid_SH','L_av_SH','L_500_SH','r_mid_SH','r_500_SH','r_av_SH',...
  'trop_SH', 'lat', 'level')
    
end



% [~,lat30] = min(abs(lat-30));
% [~,lat70] = min(abs(lat-70));
% 
% [~,p500] = min(abs(level-500));
% 
% lambda_paul = squeeze(mean(lambda_erai_NH_x(lat70:lat30,p500,:),1));
% 
% figure(1)
% plot(lambda_paul,'Linewidth',1.5); hold on
% ylabel('\lambda')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% legend('Paul'); legend boxoff
% title('Lambda 30-70 NH')