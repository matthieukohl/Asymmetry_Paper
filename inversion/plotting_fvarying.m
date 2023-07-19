clear;
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;

r=0.1;

% for i=1:16

i=3;

path_temp = strcat('/disk7/mkohl/interpolation/output1/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/u_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/v_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');


lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

event_latspan = [83:111]; %83 - 111
event_lonspan = [1:256];    % 1:256
event_timespan = [84:84];
level_indices = [5:30];
level = plevels(level_indices);
p_up = 20000;
p_w = 5000;
p_up_index = find(level==p_up);
Reduction = reduction_factor(r,level);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'temp_plevel', start, count);
v_in= ncread(path_v,'temp_plevel', start, count);
omega_in= ncread(path_omega,'temp_plevel', start, count);

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
    
    
f = 2 * Omega * sin(event_latspan / 180 * 3.1415926);
beta = 1 / R * d_dphi(f, dphi);
            
[A, B, C, J, J1, J2, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed,sigma_r,sigma_r_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, rhs_smoothed, omega_QG, omega_b, sigma_old,sigma_new,sigma_new_smoothed,del2_J] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));

for t=1:length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

% % effective static stability: sigma_eff
% temp_omega = squeeze(omega(:, :, :, t));
% [~, ~, theta_avg] = eff_stat_stab(level, T_avg(:, t), 1); % this routine is just used to calculate potential temperature                


% spatially varying static stability: sigma_accu
% this quantity is mainly used to calculate J term
Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);


%deal with points where sigma_accu is not smooth
sigma_smooth_window = 3;    
if sigma_smooth_window ~= 0
    for k = 1 : length(level)
        sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
    end
else
    sigma_accu_smoothed(:, :, :, t) = sigma_accu(:, :, :, t);
end


% term A (Q vector)

for k = 1 : length(level)

    % probably the ug and vg field need smoothing as well...
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
end
clear('temp_x', 'temp_y', 'div_temp');

% term B (planetary vorticity)

temp = d_dp(v(:, :, :, t), level);
for j = 1 : length(lat)
    B(j, :, :, t) = f(j) * beta(j) * temp(j, :, :);
end
clear('temp');

% right hand side

rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t); 

% Inversion

%omega_QG(:, :, :, t) = SIP_diagnostic1(R, f, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            %sigma(:, t), dphi, dlambda, rhs(:,:,:,t),  omega_b(:,:,:,t));
                                                            
omega_QG(:,:,:,t) = randn(length(event_latspan), length(event_lonspan), length(level),1); 

                                                            
for s=1:20
       
 for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level)
          if omega_QG(l,m,n,t)>=0
           sigma_r(l, m, n, t) = 0*sigma(n,t); 
          else 
           sigma_r(l, m, n, t) = (1-Reduction(n))*sigma(n,t);  
          end 
      end
    end
 end

J = omega_QG.*sigma_r;


for k=1:size(level)
    
    del2_J(:,:,k,t) = spherical_laplacian(J(:,:,k,t),phi,dphi,dlambda);
    
end
 
rhs_mod = rhs+del2_J;

omega_QG(:, :, :, t) = SIP_diagnostic1(R, f, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            sigma(:, t), dphi, dlambda, rhs_mod(:,:,:,t),  omega_b(:,:,:,t));                                                        
end 
                                                                                     
end 


% save(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/random/varyingf/RHS_rhs_smooth0_iter30_day',num2str(25*i),'.mat'),'omega','omega_QG','J','del2_J','rhs');
% 
% end

figure(1)
contour(omega_QG(:,:,12,1)); 


