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
[lambda_NH_theory_mid,lambda_NH_theory_500,lambda_NH_theory_av,...
lambda_NH_modal_theory_500, ...
k_NH_mid,k_NH_500,k_NH_av,k_v_NH_mid,k_v_NH_500,k_v_NH_av,L_mid_NH,L_500_NH,L_av_SH,...
S_mid_NH,S_500_NH,S_av_NH,...
r_mid_NH,r_500_NH,r_av_NH,...
r_mid_cond_NH,r_500_cond_NH,r_av_cond_NH,...
trop_NH] = deal(zeros(12,1));
[spec_NH_mid,spec_NH_500,spec_NH_av,spec_v_NH_mid,spec_v_NH_500,spec_v_NH_av] = deal(zeros(1440,12));
else
[lambda_SH_theory_mid,lambda_SH_theory_500,lambda_SH_theory_av,...
lambda_SH_modal_theory_500,...
k_SH_mid,k_SH_500,k_SH_av,k_v_SH_mid,k_v_SH_500,k_v_SH_av,L_mid_SH,L_500_SH,L_av_SH,...
S_mid_SH,S_500_SH,S_av_SH,r_mid_SH,r_500_SH,r_av_SH,trop_SH,r_mid_cond_SH,r_500_cond_SH,r_av_cond_SH] = deal(zeros(12,1));
[spec_SH_mid,spec_SH_500,spec_SH_av,spec_v_SH_mid,spec_v_SH_500,spec_v_SH_av] = deal(zeros(1440,12));
end

tic
for ii = 1:length(months)
ii
season = months{ii};

norm = 0;
lambda_mid = 0; lambda_av = 0; lambda_500 = 0; lambda_500_mode = 0;
k_mid = 0; k_av = 0; k_500 = 0;
k_v_mid = 0; k_v_av = 0; k_v_500 = 0;
spec_mid = 0; spec_av = 0; spec_500 = 0;
spec_v_mid = 0; spec_v_av = 0; spec_v_500 = 0;
r_mid = 0; r_av = 0; r_500 = 0;
r_mid_cond = 0; r_av_cond = 0; r_500_cond = 0;
S_mid = 0; S_av = 0; S_500 = 0;
L_mid = 0; L_av = 0; L_500 = 0;
trop = 0;

for year = years

disp('year:')
disp(year)

if tf_hemisphere ==1
file_v = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/v_2009_2018_NH.nc'];
file_temp = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/temp_2009_2018_NH.nc'];
file_omega = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/omega_2009_2018_NH.nc'];
else
file_v = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/v_2009_2018_SH.nc'];
file_temp = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/temp_2009_2018_SH.nc'];
file_omega = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/omega_2009_2018_SH.nc'];
end

lat = ncread(file_v,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(file_v,'longitude');lon = double(lon);

% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

time = ncread(file_v,'time');

time_range = season_selector(file_v, season, start_year, end_year);

tic

for t = time_range

% start = [1, 1, 1, t];
% count = [length(lon), length(lat), length(level), 1];

start = [1, 1, t];
count = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega_in = ncread(file_omega,'w',start,count);
T_in = ncread(file_temp,'t',start,count);
v_in = ncread(file_v,'v',start,count);

T = rearrange_era_2d(T_in,lat,lon);
omega = rearrange_era_2d(omega_in,lat,lon);
v = rearrange_era_2d(v_in,lat,lon);

clear('omega_in','T_in','v_in');


% find the 1-d zonal spectrum of omega and v

% [k_av_new,spec_av_new] = Spectrum2(omega(:,:,plow:trop_new),lat,lon);
% [k_mid_new,spec_mid_new] = Spectrum2(omega(:,:,pmid),lat,lon);
[k_500_new,spec_500_new] = Spectrum2(omega(:,:),lat,lon);

% k_av = k_av+k_av_new; spec_av = spec_av + spec_av_new;
% k_mid = k_mid+k_mid_new; spec_mid = spec_mid + spec_mid_new;
k_500 = k_500+k_500_new; spec_500 = spec_500 + spec_500_new;



% [k_v_av_new,spec_v_av_new] = Spectrum2(v(:,:,plow:trop_new),lat,lon);
% [k_mid_v_new,spec_v_mid_new] = Spectrum2(v(:,:,pmid),lat,lon);
[k_v_500_new,spec_v_500_new] = Spectrum2(v(:,:),lat,lon);

% k_v_av = k_v_av+k_v_av_new; spec_v_av = spec_v_av + specv__av_new;
% k_v_mid = k_v_mid+k_v_mid_new; spec_v_mid = spec_v_mid + spec_v_mid_new;
k_v_500 = k_v_500+k_v_500_new; spec_v_500 = spec_v_500 + spec_v_500_new;

% % output the static stability (paul's method of nondimensionalizing k)
% 
% [S_mid_new,S_500_new,S_av_new] = Static_Stability(sigma_accu,level,trop_new);
% 
% S_mid = S_mid+S_mid_new;
% 
% S_500 = S_500+S_500_new;
% 
% S_av = S_av + S_av_new;


% find the r-value to use

% [r_mid_new,r_500_new,r_av_new] = reduction_factor(r,level,trop_new);
% 
% r_mid = r_mid+r_mid_new;
% 
% r_500 = r_500+r_500_new;
% 
% r_av = r_av + r_av_new;

norm = norm + 1;


end

end



k_mid = k_mid/norm; k_500 = k_500/norm; k_av = k_av/norm;
spec_mid = squeeze(spec_mid)/norm; spec_500 = squeeze(spec_500)/norm;
spec_av = squeeze(spec_av)/norm;


k_v_mid = k_v_mid/norm; k_v_500 = k_500/norm; k_v_av = k_v_av/norm;
spec_v_mid = squeeze(spec_v_mid)/norm; spec_v_500 = squeeze(spec_v_500)/norm;
spec_v_av = squeeze(spec_v_av)/norm;


% S_mid = S_mid/norm; S_500 = S_500/norm; S_av = S_av/norm;
% 
% r_mid = r_mid/norm; r_500 = r_500/norm; r_av = r_av/norm;

% save variables

if tf_hemisphere ==1

k_NH_mid(ii) = k_mid; k_NH_500(ii) = k_500; k_NH_av(ii) = k_av;
spec_NH_mid(:,ii) = spec_mid; spec_NH_500(:,ii) = spec_500;
spec_NH_av(:,ii) = spec_av;

k_v_NH_mid(ii) = k_v_mid; k_v_NH_500(ii) = k_v_500; k_v_NH_av(ii) = k_v_av;
spec_v_NH_mid(:,ii) = spec_v_mid; spec_v_NH_500(:,ii) = spec_v_500;
spec_v_NH_av(:,ii) = spec_v_av;

% S_mid_NH(ii) = S_mid; S_500_NH(ii) = S_500; S_av_NH(ii) = S_av;
% 
% r_mid_NH(ii) = r_mid; r_500_NH(ii) = r_500; r_av_NH(ii) = r_av;
% 
else

k_SH_mid(ii) = k_mid; k_SH_500(ii) = k_500; k_SH_av(ii) = k_av;
spec_SH_mid(:,ii) = spec_mid; spec_SH_500(:,ii) = spec_500;
spec_SH_av(:,ii) = spec_av;

k_v_SH_mid(ii) = k_v_mid; k_v_SH_500(ii) = k_v_500; k_v_SH_av(ii) = k_v_av;
spec_v_SH_mid(:,ii) = spec_v_mid; spec_v_SH_500(:,ii) = spec_v_500;
spec_v_SH_av(:,ii) = spec_v_av;


% S_mid_SH(ii) = S_mid; S_500_SH(ii) = S_500; S_av_SH(ii) = S_av;
% 
% r_mid_SH(ii) = r_mid; r_500_SH(ii) = r_500; r_av_SH(ii) = r_av;

end

end

% plot the 1-D spectra of w^2 and v^2

R_eff = cosd(45)*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];

Ld = 1e6/(2*sqrt(2));

k = k*Ld; % nondimensionalize k

% define winter and summer spectra

spec_w_NH_500_winter = mean(spec_NH_500(:,[12,1,2]),2);
spec_w_NH_500_summer = mean(spec_NH_500(:,[6,7,8]),2);

spec_v_NH_500_winter = mean(spec_v_NH_500(:,[12,1,2]),2);
spec_v_NH_500_summer = mean(spec_v_NH_500(:,[6,7,8]),2);

% scaling theory for w based on omega equation 
% assume (rw)_xx - w = RHS \approx = 2 v_{xx}
% -(rk^2 + 1) spec(w) = - 2k^2 spec(v)
% spec(w)^2 = 4 k^4/(rk^2+1)^2 spec(v)^2

%r_winter = mean(r_500_NH([12,1,2]),2);
%r_summer = mean(r_500_NH([6,7,8]),2);
r_winter = 0.5;
r_summer = 0.1;

spec_w_theory_winter = 4*k.^4./(r_winter*k.^2+1).^2;
spec_w_theory_winter = spec_w_theory_winter'.*spec_v_NH_500_winter;
%spec_w_theory_winter = spec_w_theory_winter/max(abs(spec_w_theory_winter));

spec_w_theory_summer = 4*k.^4./(r_summer*k.^2+1).^2;
spec_w_theory_summer = spec_w_theory_summer'.*spec_v_NH_500_summer;
%spec_w_theory_summer = spec_w_theory_summer/max(abs(spec_w_theory_summer));


% make plots

figure

loglog(k,spec_w_NH_500_winter,'b'); hold on;
loglog(k,spec_w_NH_500_summer,'b--'); hold on;
loglog(k,spec_v_NH_500_winter,'g'); hold on;
loglog(k,spec_v_NH_500_summer,'g--'); hold on;
loglog(k,spec_w_theory_winter,'r'); hold on;
loglog(k,spec_w_theory_summer,'r--')
loglog(k,k.^(-3),'k--'); hold on;
loglog(k,k.^(-1),'k-.'); hold on;
xlabel('Wavevector k');
ylabel('spec');
legend('w winter','w summer','v winter','v summer','w theory winter','w theory summer','k^{-3}','k^{-1}','Location','SouthWest'); legend boxoff;
set(gca,'LineWidth',1.4)
set(gca,'FontSize',12);

figure

loglog(k,spec_w_NH_500_winter/max(abs(spec_w_NH_500_winter)),'b'); hold on;
loglog(k,spec_w_NH_500_summer/max(abs(spec_w_NH_500_summer)),'b--'); hold on;
loglog(k,spec_v_NH_500_winter/max(abs(spec_v_NH_500_winter)),'g'); hold on;
loglog(k,spec_v_NH_500_summer/max(abs(spec_v_NH_500_summer)),'g--'); hold on;
loglog(k,spec_w_theory_winter/max(abs(spec_w_theory_winter)),'r'); hold on;
loglog(k,spec_w_theory_summer/max(abs(spec_w_theory_summer)),'r--')
loglog(k,k.^(-3),'k--'); hold on;
loglog(k,k.^(-1),'k-.'); hold on;
xlabel('Wavevector k');
ylabel('spec');
legend('w winter','w summer','v winter','v summer','w theory winter','w theory summer','k^{-3}','k^{-1}','Location','SouthWest'); legend boxoff;
set(gca,'LineWidth',1.4)
set(gca,'FontSize',12);


% make theoretical plots

figure

f = @(x,r) 4*x./(r*x.^2+1).^2;

loglog(k,f(k,1),'b'); hold on; 
loglog(k,f(k,0.5),'r'); hold on; 
loglog(k,f(k,0.1),'g'); hold on; 
loglog(k,f(k,0.01),'c'); loglog(k,10*k,'k--');
hold on; loglog(k,1e3*k.^-3,'k-.'); 
legend('r=1','r=0.5','r=0.1','r=0.01','k','k^{-3}','Location','NorthEast'); 
legend boxoff; hold on; 
xlabel('Wavenumber k'); ylabel('Spectrum');
set(gca,'FontSize',12); 
title('\rm Theoretical Scaling formula')

f = @(x,r) 4*x.^4./(r*x.^2+1).^2;

figure

loglog(k,f(k,1).*spec_v_NH_500_winter','b'); hold on; 
loglog(k,f(k,0.3).*spec_v_NH_500_winter','r'); hold on; 
loglog(k,f(k,0.1).*spec_v_NH_500_winter','g'); hold on; 
loglog(k,f(k,0.01).*spec_v_NH_500_winter','c'); hold on;
loglog(k,10*k,'k--'); hold on; 
loglog(k,1e3*k.^-3,'k-.'); 
legend('r=1','r=0.3','r=0.1','r=0.01','k','k^{-3}','Location','NorthEast'); 
legend boxoff; hold on; 
xlabel('Wavenumber k'); ylabel('Spectrum');
set(gca,'FontSize',12); 
title('\rm Theoretical Scaling formula * spec(v^2)')




figure('Renderer', 'painters', 'Position', [10 10 900 400])

subplot(1,2,1)
loglog(k,spec_w_NH_500_winter/max(abs(spec_w_NH_500_winter)),'b'); hold on;
loglog(k,spec_v_NH_500_winter/max(abs(spec_v_NH_500_winter)),'r'); hold on;
loglog(k,k.^(-3),'k--'); hold on;
loglog(k,k.^(-1),'k-.'); hold on;
xlabel('Wavevector k');
ylabel('spec');
legend('w','v','k^{-3}','k^{-1}'); legend boxoff;
title('\rm Winter')

subplot(1,2,2)

loglog(k,spec_w_NH_500_summer/max(abs(spec_w_NH_500_summer)),'b'); hold on;
loglog(k,spec_v_NH_500_summer/max(abs(spec_v_NH_500_summer)),'r'); hold on;
loglog(k,k.^(-3),'k--'); hold on;
loglog(k,k.^(-1),'k-.'); hold on;
xlabel('Wavevector k');
ylabel('spec');
legend('w','v','k^{-3}','k^{-1}'); legend boxoff;
title('\rm Summer')







figure('Renderer', 'painters', 'Position', [10 10 900 400])

subplot(1,2,1)
loglog(k,spec_w_NH_500_winter,'b'); hold on;
loglog(k,spec_v_NH_500_winter,'r'); hold on;
loglog(k,spec_w_theory_winter,'k-.'); hold on;
loglog(k,k.^(-3),'k--'); hold on;
xlabel('Wavevector k');
ylabel('spec');
legend('w','v','theory w','k^{-3}'); legend boxoff;
title('\rm Winter')

subplot(1,2,2)

loglog(k,spec_w_NH_500_summer,'b'); hold on;
loglog(k,spec_v_NH_500_summer,'r'); hold on;
loglog(k,spec_w_theory_summer,'k-.'); hold on;
loglog(k,k.^(-3),'k--');
xlabel('Wavevector k');
ylabel('spec');
legend('w','v','theory w','k^{-3}'); legend boxoff;
title('\rm Summer')



