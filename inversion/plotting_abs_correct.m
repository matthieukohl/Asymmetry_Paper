% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (abs-simulations) iteratively using Dirichlet-BC w=0 on 
% the boundaries.

clear;

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

%list = {'abs3.0_T85','abs4.0_T85','abs6.0_T85'};

%list = {'abs1.8_T85','abs2.0_T85','abs2.5_T85'};

%list = {'abs1.2_T85','abs1.4_T85','abs1.6_T85'};

%list = {'abs0.8_T85','abs0.9_T85','abs1.0_T85'};

%list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85'};

for tt = 1:16
tt
for i= 6:6

% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/ug_day',num2str(5 *(28+i)),'0h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/vg_day',num2str(5 *(28+i)),'0h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');

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
event_timespan = [100:100];
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

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)
    
% % Smooth the T-fields
% 
% %deal with points where sigma_accu is not smooth
% sigma_smooth_window = 1;    
% if sigma_smooth_window ~= 0
%     for k = 1 : length(level)
%         T(:, :, k, t) = smooth2a(T(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
%     end
% else
%     T(:, :, :, t) = T(:, :, :, t);
% end

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

% % p_sat - saturation pressure water vapor
% 
% p_sat = p_triple*exp(-L/Rc*(1./T(:,:,:,t)-1/T_triple));
% 
% % r_sat - saturated mixing ratio
% 
% r_sat = epsilon*p_sat./(Level-p_sat); 
%  
% % Equivalent Potential Temperature (& Gradient)
% 
% theta_star = theta.*exp(L/cpc*r_sat./T(:,:,:,t));
% 
% dTheta_star_dp = d_dp(theta_star,level);
% 
% dTheta_dp = d_dp(theta,level);
% 
% % Dry-lapse rate
% 
% Gamma_dry = Ra/cp;
% 
% % Moist-lapse rate
% 
% Gamma_moist = Ra/cp*(1+L/(Ra*T(:,:,:,t)).*r_sat)./(1+(cpc/cp+(L./(Rc*T(:,:,:,t))-1)*L./(cp*T(:,:,:,t))).*r_sat);
% 
% % reduction factor r
% 
% r(:,:,:,t) = (theta./theta_star).*(Gamma_moist./Gamma_dry).*(dTheta_star_dp./dTheta_dp);
% 
% r(:,:,:,t) = max(r(:,:,:,t),0);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);

% term A (Q vector)

for k = 1 : length(level)
    
    du_dlambda(:, :, k, t) = d_dlambda(u(:, :, k, t), dlambda);
    dv_dlambda(:, :, k, t) = d_dlambda(v(:, :, k, t), dlambda);
    du_dphi(:, :, k, t) = d_dphi(u(:, :, k, t), dphi);
    dv_dphi(:, :, k, t) = d_dphi(v(:, :, k, t), dphi);

    dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
    dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
    Qx(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
        (1 ./ cos(Phi))    .* dv_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
    Qy(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                              dv_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
    div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
    A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
    
    %zeta(:, :, k, t) = curl(phi, u(:, :, k, t), v(:, :, k, t), dphi, dlambda);
end

clear('temp_x', 'temp_y', 'div_temp');

% term B (planetary vorticity)

temp = d_dp(v(:, :, :, t), level);
for j = 1 : length(lat)
    B(j, :, :, t) = f(j) * beta(j) * temp(j, :, :);
end
clear('temp');

% for j = 1 : length(lat)
%     B(j, :, :, t) = f0 * beta * temp(j, :, :);
% end
% clear('temp');

% % term C (Rayleigh-drag)
% 
% dzeta_dp = d_dp(zeta(:, :, :, t), level);
% 
% drag = 2.3*10^(-5); 
% 
% C(:,:,1,t) = -9.81*drag*dzeta_dp(:,:,1,t);

% right hand side

rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t); 

% cut to lowest non-NaN level

for bound=1:length(level)

b = isnan(rhs(:,:,bound,t)); 

if any(b(:))==0
    break
end 

end

rhs_place = rhs(:,:,bound:end,t);
omega_b_place = omega_b(:,:,bound:end,t);
sigma_place = sigma(bound:end,t);
level_place = level(bound:end);
r_place = r(:,:,bound:end,t);
omega_place = omega(:,:,bound:end,t);

% Start Iterative Inversion with Random field

%w = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
  %sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

w = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
  sigma_place(:), dphi, dlambda, rhs_place(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));

%w = SIP_diagnostic(R, f0, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
  %sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

%w = randn(length(event_latspan), length(event_lonspan), length(level_place)); 

% % Pre-Allocation
% 
% sigma_r = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));
% del2_J = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));
% 
% 
% for s=1:30
%        
%  for l=1:length(event_latspan)
%     for m=1:length(event_lonspan)
%       for n= 1:length(level_place)
%           if w(l,m,n)>=0
%            sigma_r(l, m, n) = 0*sigma_place(n); 
%           else 
%            sigma_r(l, m, n) = (1-r_place(l,m,n))*sigma_place(n);  
%           end 
%       end
%     end
%  end
% 
% J = w.*sigma_r;
% 
% 
% for k=1:size(level_place)
%     
%     del2_J(:,:,k) = spherical_laplacian(J(:,:,k),phi,dphi,dlambda);
%     
% end
%  
% rhs_mod = rhs_place+del2_J; 
% 
% %w = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%           %sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));
%                              
% %w = SIP_diagnostic(R, f0, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%                                  %sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));
%                                  
% w = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%  sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));                               
%                                                          
% end

omega_QG(:,:,bound:end,t) = w;

end

save(strcat('/net/aimsir/archive1/mkohl/inversion/output/abs_simulation/',list{tt},'/Inversion_RHS_smooth0_geostrophic_w_rand1_day',num2str(5*(28+i)),'.mat'),'omega','omega_QG','r','rhs');

end

end


figure(1)
contour(rhs(:,:,1,1)); colorbar

figure(2)
histogram(rhs(:,:,end,1),'normalization','pdf');

