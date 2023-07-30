% Compare DRV-Theory to Modal Instability Solutions

clear;
close all;

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
g = 9.80665;     % gravitational acceleration [m/s^2]
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

% Define an effective r and spatial scale b

[b_modal,R_eff,R1_eff,R2_eff,R_mean,R_500,R_300,R_700,S0_700,S0_500,S0_300,...
    S0,S0_mean,S1,S2,S0_eff,S1_eff,S2_eff,...
    H_trop,delta_p,L_D,u_eff,u_av,beta_ND,f_eff,lat_max] = ...
    deal(zeros(length(list),1));

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

for tt = 1:15
tt
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/temp_day0560h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/ug_day0560h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/vg_day0560h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/omega_day0560h00.nc');

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

event_latspan = [83:111]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [6:9];
level_indices = [1:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'ug', start, count);
v_in= ncread(path_v,'vg', start, count);
omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
    for k = 1 : length(level_indices)
        for t = 1 : length(event_timespan)
            T(:, :, k, t) = T_in(:, :, k, t)';
            u(:, :, k, t) = u_in(:, :, k, t)';
            v(:, :, k, t) = v_in(:, :, k, t)';
            omega(:, :, k, t) = omega_in(:, :, k, t)';
        end
    end
    
    
f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);
            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,lapse,density,T_virtual] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)
    
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);
density(:,:,:,t) = rho(T(:,:,:,t),Level);
T_virtual(:,:,:,t) = temp_virtual(T(:,:,:,t),Level);
lapse(:,:,:,t) = -9.81*density(:,:,:,t).*d_dp(T(:,:,:,t),level)*10^3; % K/km

% calculate the vorticity

for k = 1:length(level)
vor(:,:,k,t) = curl(phi, u(:,:,k,t)-mean(mean(u(:,:,k,t))), v(:,:,k,t)-mean(mean(v(:,:,k,t))), dphi, dlambda);
end

end

% Calculate spatial scale (b) and effective r of the dominant instability

[latm,lonm] = find(omega(:,:,12,1)==min(min(omega(:,:,12,1))));

b_modal(tt) = dlambda*R*cosd(lat(latm))*(min(abs(lonm-find(omega(latm,:,12,1)>0))));

% find troposphere

lapse_mean = mean(mean(mean(lapse,1),2),4);

trop = find(lapse_mean>-2,1);

trop_height(tt) = trop;

% find surface

surf = find(r(latm,lonm,:,1)>=0,1);

% reset r
% 
% for ii = 1:length(lat)
%     for jj = 1:length(lon)
%         for kk = 1:surf
%         for time = 1:length(event_timespan)
% if r(ii,jj,kk,time)<0
%     r(ii,jj,kk,time) = 1;
% end
%         end
%         end
%     end
% end
% 
% surf = 1;

% Calculate R_eff

dp = level(surf)-level(trop);
p = level(surf:trop);
N2 = squeeze(sigma_accu(latm,lonm,surf:trop,1));
rz = squeeze(r(latm,lonm,surf:trop,1));

%S0_300(tt) = N2(14); R_300(tt) = rz(14);
S0_500(tt) = N2(10); R_500(tt) = rz(10);
S0_700(tt) = N2(7); R_700(tt) = rz(7);

S0_eff(tt) = trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2.*rz);
S0(tt) =  trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2);
S0_mean(tt) = 1/dp*trapz(p,N2);

S1_eff(tt) = trapz(p,1/(2*dp)*(1-4/pi*sin(pi*(level(surf)-p)/dp)).*N2.*rz);
S1(tt) = trapz(p,1/(2*dp)*(1-4/pi*sin(pi*(level(surf)-p)/dp)).*N2);

S2_eff(tt) = trapz(p,1/(2*dp)*cos(pi*(level(surf)-p)/dp).*N2.*rz);
S2(tt) = trapz(p,1/(2*dp)*cos(pi*(level(surf)-p)/dp).*N2);


R_eff(tt) = S0_eff(tt)/S0(tt); R1_eff(tt) = S1_eff(tt)/S1(tt);
R2_eff(tt) = S2_eff(tt)/S2(tt); R_mean(tt) = -1/dp*trapz(p,rz);

% Tropospheric Height (cartesian coordinates)

H_trop(tt) = -Ra/g *1/dp*trapz(p,T_virtual(latm,lonm,surf:trop,1))*log(level(surf)/level(trop));

% Tropospheric Height (pressure coordinates)

delta_p(tt) = dp;

delta_p_500 = abs(level(13)-level(7));

% Effective Velocity 

u_av(tt) = abs(1/dp*trapz(p,u(latm,lonm,surf:trop,1)));

u_eff(tt) = abs(2/(sqrt(2)*dp)*trapz(p,squeeze(u(latm,lonm,surf:trop,1)).*cos(pi*(level(surf)-p)/dp)));

u_500(tt) = abs(u(13)-u(7))/2;

% Effective Coriolis

f_eff(tt) = f(latm);

% Latitude of Maximum Instability

lat_max(tt) = lat(latm);

% Deformation Radius

L_D(tt) = sqrt(-S0(tt)).*delta_p(tt)./f_eff(tt);

% Non-Dimensional beta

beta_ND(tt) = beta(latm)*L_D(tt)^2/u_eff(tt);

end

% Save Modal Variables for Theoretical Comparison

%save('/disk7/mkohl/Mode_Kohl_OGorman/GCM_Modal.mat','R_eff','R_mean','H_trop',...
    %'delta_p','surf_temp','u_eff','f_eff','S0','S0_mean','b_modal')
    
 %save('/net/halo/disk28/disk7/mkohl/Mode_Zurita_Gotor/2D_3D_data.mat',...
     %'R_eff','R_mean','R_500','R_700','beta_ND')
 
 % save('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_data.mat',...
     %'R_eff','R_mean','R_500','R_700','beta_ND')
     
%  R_alternative = R_eff;
%  
%  save('test.mat','R_alternative');
%load('test.mat')
%load('test2.mat')

%R_alternative_2 = R_eff;

%save('test2.mat','R_alternative_2')

% Calculate rescale factors DRV

% day = 24*3600; km = 1000;
% 
% rescale_sigma = u_eff.*f_eff./(delta_p.*sqrt(-S0))*day;
% 
% rescale_b = sqrt(-S0).*delta_p./f_eff*1/km;
% 
% rescale_sigma_av = u_av.*f_eff./(delta_p.*sqrt(-S0_mean))*day;
% 
% rescale_b_av = sqrt(-S0_mean).*delta_p./f_eff*1/km;

sigma_GCM = [0.2020,0.2471,0.2407,0.2667,0.2836,0.2812,0.2674,0.2299,0.2301,...
0.2147,0.2341,0.2820,0.3520,0.5490,1.1526];
% 
% sigma_DRV = [0.2453, 0.3344, 0.4003, 0.4238, 0.4400, 0.4767, 0.5035, 0.5929,...
% 0.5933,0.6472,0.6507,0.6903,0.7486,0.7492,0.8427];
% 
% sigma_DRV_av = [0.3834,0.4642,0.5208,0.5528,0.5620,0.6069,0.6443,0.7132,...
% 0.7374,0.7747,0.8014,0.8479,0.8830,0.9173,0.9923];
% 
% sigma_DRV_500 = [0.2112,0.3229,0.3808,0.4363,0.4447,0.5246,0.6043,0.7116,...
% 0.7772,0.8405,0.9182,1.0260,1.1070,1.3849,1.5476];
% 
% b_DRV = [8.7336,7.5398,6.7858,6.5345,6.3460,5.9690,5.6549,4.7124,4.7124,...
% 4.0841,4.0841,3.7071,3.0159,3.0159,2.0735];
% 
% b_DRV_av = [6.9743,6.0947,5.4664,5.0894,5.0265,4.5867,4.1469,3.3929,3.1416,...
% 2.7018,2.4504,2.0106,1.6965,1.4451,1.0053];
% 
% b_DRV_500 = [9.2363,7.6655,7.0372,6.4088,6.2832,5.4664,4.5867,3.4558,...
%  2.7018,2.0735,1.5080,0.8796,0.6283, 0.1885, 0.1257];
% 
% sigma_mode = [0.6002,0.6117,0.6262,0.6322,0.6357,0.6443,0.6505,0.6743,...
% 0.6765, 0.6978, 0.6994, 0.7165, 0.7480, 0.7484, 0.8042];
% 
% sigma_mode_av = [0.6236, 0.6417, 0.6561, 0.6651,0.6659,0.6816,0.6965,0.7296,...
% 0.7411,0.7608,0.7776,0.8064,0.8337,0.8592,0.9166];
% 
sigma_mode_beta = [0.1526,-0.0141,0.0288,0.0007,0.0016,-0.0017,-0.0067,...
0.0054,0.0034, 0.0021, -0.0009, -0.0026, 0.0022, 0.0114,-0.0004];
% 
% sigma_dry = 2-sqrt(2);
% 
% 
% b_mode = [1.6336, 1.5080,1.5080, 1.3823, 1.3823, 1.3823, 1.3823,1.1310,...
% 1.2566, 1.1310,1.0053, 1.0053, 0.8796, 0.8796, 0.7540];
% 
% b_mode_av = [1.5080,1.3823,1.3823,1.2566,1.2566,1.1310,1.1310,1.0053,1.0053,...
% 0.8796,0.8796,0.7540,0.6283,0.6283,0.5027];

% Calculate rescale factors modes

% rescale_sigma_mode = u_eff.*f_eff./(delta_p.*sqrt(-S0))*day*1/sqrt(2);
% 
% rescale_b_mode = sqrt(-S0).*delta_p./f_eff*1/km*sqrt(2);
% 
% rescale_sigma_mode_av = u_eff.*f_eff./(delta_p.*sqrt(-S0_mean))*day*1/(sqrt(2));
% 
% rescale_b_mode_av = sqrt(-S0_mean).*delta_p./f_eff*1/km*sqrt(2);

%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/DRV.mat','sigma_DRV',...
    %'sigma_DRV_av','b_DRV','b_DRV_av')
    
%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/DRV_500.mat','sigma_DRV_500',...
%'b_DRV_500')
    
    
%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode.mat','sigma_mode',...
    %'sigma_mode_av','b_mode','b_mode_av')
    
%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_500.mat','sigma_mode_500')

%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_500_beta.mat','sigma_mode_500_beta')

%save('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_beta.mat','sigma_mode_beta')

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/DRV.mat');

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/DRV_500.mat');

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode.mat');

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_500.mat');

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_500_beta.mat');

load('/disk7/mkohl/Mode_Kohl_OGorman/2D_3D_comparison/Mode_beta.mat');



[rescale_sigma,rescale_b,rescale_sigma_av,rescale_b_av,...
 rescale_sigma_500,rescale_sigma_b_500] = ...
 rescale_2D_3D(u_eff,u_av,u_500,f_eff,delta_p,delta_p_500,S0,S0_mean,S0_500);

sigma_dry = sqrt(2)-1;

    

figure(1)
plot(surf_temp,R_eff,'linewidth',1.5); hold on;
%plot(surf_temp,R_mean,'linewidth',1.5); hold on;
plot(surf_temp,R_alternative,'linewidth',1.5); hold on;
plot(surf_temp,R_alternative_2,'linewidth',1.5)
yticks([0 0.2 0.4 0.6 0.8 1.0])
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'})
xlabel('T_g (K)')
ylabel('r')
%legend('r_{eff}','r_{av}','Location','southwest')
legend('r_{b0}','r_{b1}','r_{nob}','Location','southwest')
legend boxoff
%title('Reduction Factors Mode')
title('Effective Reduction Factors Mode')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(3)
plot(surf_temp,-S0,'linewidth',1.5); hold on;
%plot(surf_temp,-S0_mean,'linewidth',1.5); hold on;
plot(surf_temp,S0_500/2,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('S')
%legend('S_0','S_{av}')
legend('S_0','S_500')
legend boxoff
title('Static Stability Mode')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

%

% 
figure(6)
plot(surf_temp,u_eff,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('u (ms)')
title('Effective Zonal Velocity')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);



figure(10)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_DRV,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode'.*sigma_mode,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode'.*sigma_dry,'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (effective r)')
legend('GCM','DRV','Mode','Dry')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(11)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_av'.*sigma_DRV_av,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode_av'.*sigma_mode_av,'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (Mean r)')
legend('GCM','DRV','Mode')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(12)
plot(surf_temp,b_modal/km,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b'.*b_DRV,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mode'.*b_mode,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('b (km)')
title('Half-Ascent Area (Effective)')
legend('GCM','DRV','Mode')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(13)
plot(surf_temp,b_modal/km,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_av'.*b_DRV_av,'linewidth',1.5); hold on;
plot(surf_temp,rescale_b_mode_av'.*b_mode_av,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('b (km)')
title('Half-Ascent Area (Mean)')
legend('GCM','DRV','Mode')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(14)
plot(surf_temp,sigma_DRV,'linewidth',1.5); hold on;
plot(surf_temp,sigma_mode/sqrt(2),'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\sigma')
title('Growthrate (ND)')
legend('DRV','Mode')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(15)
plot(surf_temp,beta_ND,'linewidth',1.5); hold on;
xlabel('T_g (K)')
ylabel('\beta')
title('Beta (ND)')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(16)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_DRV,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_mode/sqrt(2),'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_dry/sqrt(2),'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_mode_beta,'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (effective r)')
legend('GCM','DRV','Mode','Dry Mode','Mode beta')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(17)
%plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma'.*sigma_DRV,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode'.*sigma_mode,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode'.*sigma_dry,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_mode'.*sigma_mode_beta*sqrt(2),'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (effective r)')
legend('DRV','Mode','Dry Mode','Mode beta')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(20)
plot(surf_temp,R_eff,'linewidth',1.5); hold on;
plot(surf_temp,R_mean,'linewidth',1.5); hold on;
plot(surf_temp,R_500,'linewidth',1.5); hold on;
plot(surf_temp,R_700,'linewidth',1.5)
yticks([0 0.2 0.4 0.6 0.8 1.0])
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'})
xlabel('T_g (K)')
ylabel('r')
legend('r_{eff}','r_{av}','r_{500}','r_{700}','Location','southwest')
legend boxoff
%title('Reduction Factors Mode')
title('Reduction Factors Mode')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(21)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_500.*sigma_DRV_500,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_500.*sigma_mode_500,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma_500.*sigma_dry,'linewidth',1.5); hold on;
%plot(surf_temp,rescale_sigma_500.*sigma_mode_500_beta,'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (r and S at 500hPA, u_{eff}, H_{trop})')
legend('GCM','DRV','Mode','Dry Mode','Location','NorthWest')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(22)
plot(surf_temp,sigma_GCM,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma.*sigma_DRV_500,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma.*sigma_mode_500,'linewidth',1.5); hold on;
plot(surf_temp,rescale_sigma.*sigma_dry,'linewidth',1.5); hold on;
%plot(surf_temp,rescale_sigma.*sigma_mode_500_beta,'linewidth',1.5)
xlabel('T_g (K)')
ylabel('\sigma (day^{-1})')
title('Growthrate (r at 500hPA, S_{eff}, u_{eff}, H_{trop})')
legend('GCM','DRV','Mode','Dry Mode','Location','NorthWest') 
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

figure(23)
plot(squeeze(-omega(latm,lonm,:,1))/max(-omega(latm,lonm,:,1)),level/100,'linewidth',1.5); hold on;
plot(squeeze(r(latm,lonm,:,1)),level/100,'linewidth',1.5); hold on;
plot(sin(pi*(level(surf)-p)/dp),p/100,'linewidth',1.5); hold on;
plot(squeeze(sigma_accu(latm,lonm,:,1))/max(sigma_accu(latm,lonm,:,1)),level/100,'linewidth',1.5); hold on;
xlabel('Amplitude')
ylabel('p(hPA)')
legend('\omega/\omega_{max}','r','sin','N2','Location','NorthWest') 
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);
ylim([level(end)/100 level(1)/100])
set(gca,'Ydir','reverse')















