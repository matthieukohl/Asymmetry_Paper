% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (abs-simulations) iteratively using Dirichlet-BC w=0 on 
% the boundaries.

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
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

index = 11*ones(size(surf_temp));

[Norm_RHS,Norm_Deformation,Norm_A1,Norm_A1_bb,Norm_A1_pp,Norm_A1_bp,Norm_A1_pb,...
Norm_A3,Norm_zeta,Norm_zeta_bar,Norm_zeta_prime,...
Norm_zeta_x,Norm_zeta_x_bar,Norm_zeta_x_prime,...
Norm_zeta_y,Norm_zeta_y_bar,Norm_zeta_y_prime,...
Norm_zeta_px,Norm_zeta_px_bar,Norm_zeta_px_prime,...
Norm_zeta_py,Norm_zeta_py_bar,Norm_zeta_py_prime,...
Norm_du_dp,Norm_dv_dp,Norm_du_dp_bar,...
Norm_dv_dp_bar,Norm_du_dp_prime,Norm_dv_dp_prime,...
Skew_RHS,Skew_A1,Skew_A1_bb,Skew_A1_pp,Skew_A1_pb,Skew_A1_bp,...
Skew_shear_u,Skew_shear_u_bar,Skew_shear_u_prime,Skew_shear_v,Skew_shear_v_bar,...
Skew_shear_v_prime,Skew_vor,Skew_vor_bar,Skew_vor_prime,...
Skew_vor_x,Skew_vor_x_bar,Skew_vor_x_prime,...
Skew_vor_y,Skew_vor_y_bar,Skew_vor_y_prime,...
Skew_vor_px,Skew_vor_px_bar,Skew_vor_px_prime,...
Skew_vor_py,Skew_vor_py_bar,Skew_vor_py_prime,...
Skew_Deformation,Skew_Product,Skew_Forcing,Skew_Forcing_Product] = deal(zeros(length(surf_temp),1));

for tt = 14:14
    
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
event_timespan = [6:9]; % 6-9
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

[T, u, v, omega,vor,zeta,zeta_x,zeta_y,zeta_px,zeta_py,Forcing] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
    for k = 1 : length(level_indices)
        for t = 1 : length(event_timespan)
            T(:, :, k, t) = T_in(:, :, k, t)';
            u(:, :, k, t) = u_in(:, :, k, t)';
            v(:, :, k, t) = v_in(:, :, k, t)';
            omega(:, :, k, t) = omega_in(:, :, k, t)';
        end
    end
    
    
f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, J, Qx, Qy, ...
 du_dlambda, dv_dlambda, du_dphi, dv_dphi,du_dp,dv_dp,zeta_p, ...
 dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)
    
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);



for k = 1:length(level)
    
    zeta(:,:,k,t) = curl(phi,u(:,:,k,t),v(:,:,k,t),dphi,dlambda);
    
end

du_dp(:,:,:,t) = d_dp(u(:,:,:,t),level);
dv_dp(:,:,:,t) = d_dp(v(:,:,:,t),level);
zeta_p(:,:,:,t) = d_dp(zeta(:, :, :, t), level);

for k = 1:length(level)
    
zeta_x(:,:,k,t) = d_dx(phi, zeta(:,:,k,t), dlambda);
zeta_y(:,:,k,t) = d_dy(phi, zeta(:,:,k,t), dphi);
zeta_px(:,:,k,t) = d_dx(phi, zeta_p(:,:,k,t), dlambda);
zeta_py(:,:,k,t) = d_dy(phi, zeta_p(:,:,k,t), dlambda);
Forcing(:,:,k,t) = v_del(phi, du_dp(:, :, k,t), dv_dp(:, :, k,t), zeta_p(:, :, k, t), dphi, dlambda);
    
end

Forcing_Product = du_dp.*zeta_px;

% for k = 1:length(level)
%     
%     zeta_x(:,:,k,t) = d_dx(phi, zeta(:,:,k,t), dlambda);
%     zeta_y(:,:,k,t) = d_dy(phi, zeta(:,:,k,t), dphi);
%     zeta_px(:,:,k,t) = d_dx(phi, zeta_p(:,:,k,t), dlambda);
%     zeta_py(:,:,k,t) = d_dx(phi, zeta_p(:,:,k,t), dlambda);
%     
% end

rr = r_parameter1(T(:,:,:,t),Level,level);
rr(omega(:,:,:,t)>0) = 1;
%diab(:,:,:,t) = (1-rr).*omega(:,:,:,t).*S(:,:,:,t);
%diab(:,:,:,t) = diab(:,:,:,t).*(1e5./Level).^kappa;


end

[A, A1, A2, A3] = traditional_A(phi,level,u,v,zeta,event_timespan,dphi,dlambda,f0,beta);

B = traditional_B(phi,level,u,v,T,event_timespan,dphi,dlambda);



Deform = Deformation(phi,level, u, v, T, event_timespan, dphi, dlambda); %A+B - 2*A1-A3;

rhs = 2*A1+A3+Deform;

u_prime = u -mean(u(:,:,:,:),2);
v_prime = v-mean(v(:,:,:,:),2);

u_bar = u-u_prime;
v_bar = v-v_prime;

zeta_prime = zeta - mean(zeta(:,:,:,:),2); zeta_bar = zeta - zeta_prime;

zeta_x_prime = zeta_x - mean(zeta_x(:,:,:,:),2); zeta_x_bar = zeta_x - zeta_x_prime;

zeta_y_prime = zeta_y - mean(zeta_y(:,:,:,:),2); zeta_y_bar = zeta_y - zeta_y_prime;

zeta_px_prime = zeta_px - mean(zeta_px(:,:,:,:),2); zeta_px_bar = zeta_px - zeta_px_prime;

zeta_py_prime = zeta_py - mean(zeta_py(:,:,:,:),2); zeta_py_bar = zeta_py - zeta_py_prime;

du_dp_prime = du_dp-mean(du_dp(:,:,:,:),2);

dv_dp_prime = dv_dp - mean(dv_dp(:,:,:,:),2);

du_dp_bar = du_dp - du_dp_prime;

dv_dp_bar = dv_dp - dv_dp_prime;

Product = du_dp_bar.*zeta_x_prime;

[A_bb, A1_bb, A2_bb, A3_bb] = traditional_A(phi,level,u_bar,v_bar,zeta_bar,event_timespan,dphi,dlambda,f0,beta);

[A_pp, A1_pp, A2_pp, A3_pp] = traditional_A(phi,level,u_prime,v_prime,zeta_prime,event_timespan,dphi,dlambda,f0,beta);

[A_bp, A1_bp, A2_bp, A3_bp] = traditional_A(phi,level,u_bar,v_bar,zeta_prime,event_timespan,dphi,dlambda,f0,beta);

[A_pb, A1_pb, A2_pb, A3_pb] = traditional_A(phi,level,u_prime,v_prime,zeta_bar,event_timespan,dphi,dlambda,f0,beta);

%Norm of Different Terms 500hPA

[Norm_RHS_p,Norm_Deformation_p,Norm_A1_p,Norm_A1_bb_p,Norm_A1_pp_p,...
 Norm_A1_bp_p,Norm_A1_pb_p,Norm_A3_p,Norm_zeta_p,Norm_zeta_bar_p,Norm_zeta_prime_p,...
 Norm_du_dp_p,Norm_du_dp_bar_p,Norm_du_dp_prime_p,Norm_dv_dp_p,Norm_dv_dp_bar_p,...
 Norm_dv_dp_prime_p,Norm_zeta_x_p,Norm_zeta_x_bar_p,Norm_zeta_x_prime_p,...
 Norm_zeta_y_p,Norm_zeta_y_bar_p,Norm_zeta_y_prime_p,...
 Norm_zeta_px_p,Norm_zeta_px_bar_p,Norm_zeta_px_prime_p,...
 Norm_zeta_py_p,Norm_zeta_py_bar_p,Norm_zeta_py_prime_p,...
 ] = deal(zeros(length(event_timespan),1));

for t = 1:length(event_timespan)

rhs_term = rhs(:,:,index(tt),t); Norm_RHS_p(t) = norm(rhs_term(:));    
    
Deformation_term(:,:) = Deform(:,:,index(tt),t); Norm_Deformation_p(t) = norm(Deformation_term(:));

A1_term = 2*A1(:,:,index(tt),t); Norm_A1_p(t) = norm(A1_term(:));

A1_bb_term = 2*A1_bb(:,:,index(tt),t); Norm_A1_bb_p(t) = norm(A1_bb_term(:));

A1_pp_term = 2*A1_pp(:,:,index(tt),t); Norm_A1_pp_p(t) = norm(A1_pp_term(:));

A1_bp_term = 2*A1_bp(:,:,index(tt),t); Norm_A1_bp_p(t) = norm(A1_bp_term(:));

A1_pb_term = 2*A1_pb(:,:,index(tt),t); Norm_A1_pb_p(t) = norm(A1_pb_term(:));


A3_term = A3(:,:,index(tt),t); Norm_A3_p(t) = norm(A3_term(:));

zeta_term = zeta(:,:,index(tt),t); Norm_zeta_p(t) = norm(zeta_term(:));

zeta_bar_term = zeta_bar(:,:,index(tt),t); Norm_zeta_bar_p(t) = norm(zeta_bar_term(:));

zeta_prime_term = zeta_prime(:,:,index(tt),t); Norm_zeta_prime_p(t) = norm(zeta_prime_term(:));

zeta_x_bar_term = zeta_x_bar(:,:,index(tt),t); Norm_zeta_x_bar_p(t) = norm(zeta_x_bar_term(:));

zeta_x_prime_term = zeta_x_prime(:,:,index(tt),t); Norm_zeta_x_prime_p(t) = norm(zeta_x_prime_term(:));

zeta_y_bar_term = zeta_y_bar(:,:,index(tt),t); Norm_zeta_y_bar_p(t) = norm(zeta_y_bar_term(:));

zeta_y_prime_term = zeta_y_prime(:,:,index(tt),t); Norm_zeta_y_prime_p(t) = norm(zeta_y_prime_term(:));

zeta_px_bar_term = zeta_px_bar(:,:,index(tt),t); Norm_zeta_px_bar_p(t) = norm(zeta_px_bar_term(:));

zeta_px_prime_term = zeta_px_prime(:,:,index(tt),t); Norm_zeta_px_prime_p(t) = norm(zeta_px_prime_term(:));

zeta_py_bar_term = zeta_py_bar(:,:,index(tt),t); Norm_zeta_py_bar_p(t) = norm(zeta_py_bar_term(:));

zeta_py_prime_term = zeta_py_prime(:,:,index(tt),t); Norm_zeta_py_prime_p(t) = norm(zeta_py_prime_term(:));

du_dp_term = du_dp(:,:,index(tt),t); Norm_du_dp_p(t) = norm(du_dp_term(:));

du_dp_bar_term = du_dp_bar(:,:,index(tt),t); Norm_du_dp_bar_p(t) = norm(du_dp_bar_term(:));

du_dp_prime_term = du_dp_prime(:,:,index(tt),t); Norm_du_dp_prime_p(t) = norm(du_dp_prime_term(:));

dv_dp_term = dv_dp(:,:,index(tt),t); Norm_dv_dp_p(t) = norm(dv_dp_term(:));

dv_dp_bar_term = dv_dp_bar(:,:,index(tt),t); Norm_dv_dp_bar_p(t) = norm(dv_dp_bar_term(:));

dv_dp_prime_term = dv_dp_prime(:,:,index(tt),t); Norm_dv_dp_prime_p(t) = norm(dv_dp_prime_term(:));

end

Norm_RHS(tt) = mean(Norm_RHS_p);

Norm_Deformation(tt) = mean(Norm_Deformation_p);

Norm_A1(tt) = mean(Norm_A1_p);

Norm_A1_bb(tt) = mean(Norm_A1_bb_p); Norm_A1_pp(tt) = mean(Norm_A1_pp_p);

Norm_A1_bp(tt) = mean(Norm_A1_bp_p); Norm_A1_pb(tt) = mean(Norm_A1_pb_p);

Norm_A3(tt) = mean(Norm_A3_p);

Norm_zeta(tt) = mean(Norm_zeta_p); Norm_zeta_bar(tt) = mean(Norm_zeta_bar_p); Norm_zeta_prime(tt) = mean(Norm_zeta_prime_p);

Norm_zeta_x(tt) = mean(Norm_zeta_x_p); Norm_zeta_x_bar(tt) = mean(Norm_zeta_x_bar_p); Norm_zeta_x_prime(tt) = mean(Norm_zeta_x_prime_p);

Norm_zeta_y(tt) = mean(Norm_zeta_y_p); Norm_zeta_y_bar(tt) = mean(Norm_zeta_y_bar_p); Norm_zeta_y_prime(tt) = mean(Norm_zeta_y_prime_p);

Norm_zeta_px(tt) = mean(Norm_zeta_px_p); Norm_zeta_px_bar(tt) = mean(Norm_zeta_px_bar_p); Norm_zeta_px_prime(tt) = mean(Norm_zeta_px_prime_p);

Norm_zeta_py(tt) = mean(Norm_zeta_py_p); Norm_zeta_py_bar(tt) = mean(Norm_zeta_py_bar_p); Norm_zeta_py_prime(tt) = mean(Norm_zeta_py_prime_p);

Norm_du_dp(tt) = mean(Norm_du_dp_p); Norm_du_dp_bar(tt) = mean(Norm_du_dp_bar_p); Norm_du_dp_prime(tt) = mean(Norm_du_dp_prime_p);

Norm_dv_dp(tt) = mean(Norm_dv_dp_p); Norm_dv_dp_bar(tt) = mean(Norm_dv_dp_bar_p); Norm_dv_dp_prime(tt) = mean(Norm_dv_dp_prime_p);

%Skewness of Different Terms 500hPA

Skew_RHS(tt) =  -Skewness(rhs(:,:,index(tt),:));

%Skew_diab(tt) = -Skewness(diab(:,:,index(tt),:));

Skew_A1(tt) = -Skewness(A1(:,:,index(tt),:));

Skew_A1_pp(tt) = -Skewness(A1_pp(:,:,index(tt),:)); Skew_A1_bb(tt) = -Skewness(A1_bb(:,:,index(tt),:));

Skew_A1_pb(tt) = -Skewness(A1_pb(:,:,index(tt),:)); Skew_A1_bp(tt) = -Skewness(A1_bp(:,:,index(tt),:));

Skew_shear_u(tt) = -Skewness(du_dp(:,:,index(tt),:)); 

Skew_shear_u_bar(tt) = -Skewness(du_dp_bar(:,:,index(tt),:)); Skew_shear_u_prime(tt) = -Skewness(du_dp_prime(:,:,index(tt),:));

Skew_shear_v(tt) = -Skewness(dv_dp(:,:,index(tt),:));

Skew_shear_v_bar(tt) = -Skewness(dv_dp_bar(:,:,index(tt),:)); Skew_shear_v_prime(tt) = -Skewness(dv_dp_prime(:,:,index(tt),:));

Skew_vor(tt) = -Skewness(zeta(:,:,index(tt),:)); Skew_vor_bar(tt) = -Skewness(zeta_bar(:,:,index(tt),:));

Skew_vor_prime(tt) = -Skewness(zeta_prime(:,:,index(tt),:));

Skew_vor_x(tt) = -Skewness(zeta_x(:,:,index(tt),:)); Skew_vor_x_bar(tt) = -Skewness(zeta_x_bar(:,:,index(tt),:)); Skew_vor_x_prime(tt) = -Skewness(zeta_x_prime(:,:,index(tt),:));

Skew_vor_y(tt) = -Skewness(zeta_y(:,:,index(tt),:)); Skew_vor_y_bar(tt) = -Skewness(zeta_y_bar(:,:,index(tt),:)); Skew_vor_y_prime(tt) = -Skewness(zeta_y_prime(:,:,index(tt),:));

Skew_vor_px(tt) = -Skewness(zeta_px(:,:,index(tt),:)); Skew_vor_px_bar(tt) = -Skewness(zeta_px_bar(:,:,index(tt),:)); Skew_vor_px_prime(tt) = -Skewness(zeta_px_prime(:,:,index(tt),:));

Skew_vor_py(tt) = -Skewness(zeta_py(:,:,index(tt),:)); Skew_vor_py_bar(tt) = -Skewness(zeta_py_bar(:,:,index(tt),:)); Skew_vor_py_prime(tt) = -Skewness(zeta_py_prime(:,:,index(tt),:));



Skew_Deformation(tt) = -Skewness(Deform(:,:,index(tt),:));

Skew_Product(tt) = -Skewness(Product(:,:,index(tt),:));

Skew_Forcing(tt) = -Skewness(Forcing(:,:,index(tt),:));

Skew_Forcing_Product(tt) = -Skewness(Forcing_Product(:,:,index(tt),:));



%Lambda_parameter(tt) = Lambda(A1(:,:,index(tt),:));

end
close all;
%figure('pos', [10, 10, 500, 300])
%set(gcf,'DefaultAxesFontSize', 12)
figure(1)
plot(surf_temp,Skew_RHS); hold on;
plot(surf_temp,Skew_A1); hold on;
plot(surf_temp,Skew_Deformation)
xlabel('T_{g}')
ylabel('Skew')
title('Skew RHS vs. 2*J(\tau,del \phi) vs. Deformation at 500hPA for abs-mode')
legend('RHS','Vorticity Advection','Deformation','location','northwest')
legend boxoff 
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
%savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_1.fig')
close all;
figure(2)
%set(gcf,'DefaultAxesFontSize', 14)
plot(surf_temp,Skew_A1); hold on;
plot(surf_temp,Skew_shear_u); hold on;
plot(surf_temp,Skew_shear_v); hold on;
plot(surf_temp,Skew_vor); hold on;
plot(surf_temp,Skew_vor_x); hold on;
plot(surf_temp,Skew_vor_y); hold on;
plot(surf_temp,Skew_Product); hold on;
xlabel('T_{g}')
ylabel('Skew')
title('Skew J(\tau,del \phi) Decomposition at 500hPA for abs-mode')
legend('Total','du/dp','dv/dp','Zeta','Zeta_x','Zeta_y','Product','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_2.fig')
% close all;
figure(3)
plot(surf_temp,Skew_A1); hold on;
plot(surf_temp,Skew_A1_bb); hold on;
plot(surf_temp,Skew_A1_pp); hold on;
plot(surf_temp,Skew_A1_pb); hold on;
plot(surf_temp,Skew_A1_bp)
xlabel('T_{g}')
ylabel('Skew')
title('Skew J(\tau,del \phi) Decomposition at 500hPA for abs-mode')
legend('A1','A1_{bb}','A1_{pp}','A1_{pb}','A1_{bp}','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_3.fig')
% close all;
figure(10)
plot(surf_temp,Skew_shear_u); hold on;
plot(surf_temp,Skew_shear_u_bar); hold on;
plot(surf_temp,Skew_shear_u_prime); hold on;
xlabel('T_{g}')
ylabel('Skew')
title('Skew du/dp Decomposition at 500hPA for abs-mode')
%legend('Full','Prime','location','northwest')
legend('Full','Bar','Prime','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_4.fig')
% close all;
% figure
% plot(surf_temp,Skew_shear_v); hold on;
% plot(surf_temp,Skew_shear_v_bar); hold on;
% plot(surf_temp,Skew_shear_v_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew dv/dp Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_5.fig')
% close all;
% figure
% plot(surf_temp,Skew_vor); hold on;
% plot(surf_temp,Skew_vor_bar); hold on;
% plot(surf_temp,Skew_vor_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew Zeta Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_6.fig')
% close all;
% figure
% plot(surf_temp,Skew_vor_x); hold on;
% plot(surf_temp,Skew_vor_x_bar); hold on;
% plot(surf_temp,Skew_vor_x_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew Zeta_x Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_7.fig')
% close all;
% figure
% plot(surf_temp,Skew_vor_y); hold on;
% plot(surf_temp,Skew_vor_y_bar); hold on;
% plot(surf_temp,Skew_vor_y_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew Zeta_y Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_8.fig')
% close all;
% 
% figure
% plot(surf_temp,Skew_Forcing); hold on;
% plot(surf_temp,Skew_Forcing_Product); hold on;
% plot(surf_temp,Skew_shear_u); hold on;
% plot(surf_temp,Skew_shear_v); hold on;
% plot(surf_temp,Skew_vor_px); hold on;
% plot(surf_temp,Skew_vor_py); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew J(\tau,del \tau) Decomposition at 500hPA for abs-mode')
% legend('Forcing','Product','du/dp','dv/dp','zeta_{px}','zeta_{py}','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_9.fig')
% close all;
% 
% figure
% plot(surf_temp,Skew_vor_px); hold on;
% plot(surf_temp,Skew_vor_px_bar); hold on;
% plot(surf_temp,Skew_vor_px_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew zeta_{px} Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_10.fig')
% close all;
% 
% figure
% plot(surf_temp,Skew_vor_py); hold on;
% plot(surf_temp,Skew_vor_py_bar); hold on;
% plot(surf_temp,Skew_vor_py_prime); hold on;
% xlabel('T_{g}')
% ylabel('Skew')
% title('Skew zeta_{py} Decomposition at 500hPA for abs-mode')
% legend('Full','Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_mode_11.fig')
% close all;
% 
% 
figure(4)
%set(gcf,'DefaultAxesFontSize', 14)
semilogy(surf_temp,Norm_RHS); hold on;
semilogy(surf_temp,Norm_A1); hold on;
semilogy(surf_temp,Norm_A3); hold on;
semilogy(surf_temp,Norm_Deformation); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('RHS-terms at 500hPA for abs-mode')
legend('RHS','Vort. Advection','Beta','Deformation','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_1.fig')
% close all;
% figure(5)
% %set(gcf,'DefaultAxesFontSize', 12)
% %semilogy(surf_temp,Norm_zeta); hold on;
% semilogy(surf_temp,Norm_zeta_bar); hold on;
% semilogy(surf_temp,Norm_zeta_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Zeta Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_2.fig')
% close all;
% figure(6)
% %set(gcf,'DefaultAxesFontSize', 14)
% %semilogy(surf_temp,Norm_du_dp); hold on;
% semilogy(surf_temp,Norm_du_dp_bar); hold on;
% semilogy(surf_temp,Norm_du_dp_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Shear u-Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northeast')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_3.fig')
% close all;
% figure(7)
% %set(gcf,'DefaultAxesFontSize', 14)
% %semilogy(surf_temp,Norm_du_dp); hold on;
% semilogy(surf_temp,Norm_dv_dp_bar); hold on;
% semilogy(surf_temp,Norm_dv_dp_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Shear v-Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northeast')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_4.fig')
% close all;
figure(8)
%set(gcf,'DefaultAxesFontSize', 14)
semilogy(surf_temp,Norm_A1); hold on;
semilogy(surf_temp,Norm_A1_bb); hold on;
semilogy(surf_temp,Norm_A1_pp); hold on;
semilogy(surf_temp,Norm_A1_pb); hold on;
semilogy(surf_temp,Norm_A1_bp); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('A Decomposition at 500hPA for abs-mode')
legend('Total','A1_{bb}','A1_{pp}','A1_{pb}','A1_{bp}','location','northeast')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_5.fig')
% 
% figure(9)
% %set(gcf,'DefaultAxesFontSize', 12)
% %semilogy(surf_temp,Norm_zeta); hold on;
% semilogy(surf_temp,Norm_zeta_x_bar); hold on;
% semilogy(surf_temp,Norm_zeta_x_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Zeta_x Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_6.fig')
% close all;
% 
% 
% figure(10)
% %set(gcf,'DefaultAxesFontSize', 12)
% %semilogy(surf_temp,Norm_zeta); hold on;
% semilogy(surf_temp,Norm_zeta_y_bar); hold on;
% semilogy(surf_temp,Norm_zeta_y_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Zeta_y Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_7.fig')
% close all;
% 
% figure(10)
% %set(gcf,'DefaultAxesFontSize', 12)
% %semilogy(surf_temp,Norm_zeta); hold on;
% semilogy(surf_temp,Norm_zeta_px_bar); hold on;
% semilogy(surf_temp,Norm_zeta_px_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Zeta_{px} Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_8.fig')
% close all;
% 
% figure(11)
% %set(gcf,'DefaultAxesFontSize', 12)
% %semilogy(surf_temp,Norm_zeta); hold on;
% semilogy(surf_temp,Norm_zeta_py_bar); hold on;
% semilogy(surf_temp,Norm_zeta_py_prime); hold on;
% xlabel('T_{g}')
% ylabel('Norm')
% title('Zeta_{py} Decomposition at 500hPA for abs-mode')
% legend('Bar','Prime','location','northwest')
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% legend boxoff 
% set(gca,'linewidth',1.5)
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_mode_9.fig')
% close all;
% 
