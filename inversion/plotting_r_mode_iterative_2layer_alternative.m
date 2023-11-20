% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (r-simulations) iteratively using Dirichlet-BC w=0 
% on the boundaries.

clear;

list = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};
 
R_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];
tic
for tt = 1:10
    
tt    

% Define Reduction Factor

r=R_factor(tt); 

% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% Define Path 

path_temp = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','temp_day0560h00.nc');
path_u = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','ug_day0560h00.nc');
path_v = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','vg_day0560h00.nc');
path_omega = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','omega_day0560h00.nc');

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
event_timespan = [191:194]; %191:194
level_indices = [1:30]; % 5-30;
%level_indices = [5,10,17];
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
    
% varying f    
% f = 2 * Omega * sind(lat);
% beta = 1 / R * d_dphi(f, dphi);

% constant f0
lat0 = lat(1)+(lat(end)-lat(1))/2;
f = 2 * Omega * sind(lat0*ones(size(lat)));
beta = 2 * Omega * cosd(lat0*ones(size(lat)))/R;

            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), 3, length(event_timespan)));

err_rms  = zeros(length(event_timespan),1);
for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

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

end

% extract only few layers

indx = [5,10,17];
rhs = rhs(:,:,indx,:);
omega = omega(:,:,indx,:);
omega_b  = omega_b(:,:,indx,:);
%r = r(:,:,indx,:);
level = level(indx);
sigma = sigma(indx,:);

for t = 1:length(event_timespan)

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
%omega_place = omega(:,:,bound:end,t);
omega_place = omega_b(:,:,bound:end,t);

% Start Iterative Inversion with w_QG dry or Random field


w_old = randn(length(lat),length(lon),length(level_place));

w_new = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
        sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

% w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
% sigma_place(:), dphi, dlambda, rhs_place(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));
    
diff = w_new-w_old;

% Pre-Allocation

sigma_r = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));
del2_J = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));

count = 0;

% Define Reduction Factor

Reduction = reduction_factor(r,level_place);

while rms(diff(:))>=10^(-6) && count<400

count = count + 1;
       
 for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level_place)
          if w_new(l,m,n)>=0
           sigma_r(l, m, n) = 0*sigma_place(n); 
          else 
           sigma_r(l, m, n) = (1-Reduction(n))*sigma_place(n);  
          end 
      end
    end
 end

J = w_new.*sigma_r;


for k=1:size(level_place)
    
    del2_J(:,:,k) = spherical_laplacian(J(:,:,k),phi,dphi,dlambda);
    
end
 
rhs_mod = rhs_place+del2_J; 

w_old = w_new;

w_new = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
           sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));
       
% w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%   sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));

diff = w_new-w_old;  

error(count) = rms(diff(:));
asym(count) = Lambda(w_new(:,:,2));

end

omega_QG(:,:,bound:end,t) = w_new;
err_rms(t) = rms(diff(:));

% figure('Renderer', 'painters', 'Position', [10 10 1200 600])
% 
% subplot(2,1,1)
% %cmax = 1; 
% cmax = max(max(-omega_QG(:,:,2,t)));
% nint = 20;
% cint = cmax/nint;
% 
% contourf(-omega_QG(:,:,2,t),-cmax:cint:cmax,'EdgeColor','none');
% colorbar
% colormap(redblue(nint));
% caxis([-cmax cmax])
% 
% subplot(2,1,2)
% 
% %cmax = 1; 
% cmax = max(max(-omega(:,:,2,t)));
% nint = 20;
% cint = cmax/nint;
% 
% contourf(-omega(:,:,2,t),-cmax:cint:cmax,'EdgeColor','none');
% colorbar
% colormap(redblue(nint));
% caxis([-cmax cmax])


end

%save(strcat('/net/aimsir/archive1/mkohl/inversion/output/r_simulation_mode/rescale',list{tt},'/final/LHS_rhs_smooth0_rms6_geostrophic_day0560h00.mat'),'omega','omega_QG','rhs');
save(strcat('/net/graupel/archive1/users/mkohl/inversion/output/r_mode_simulation/rescale',list{tt},'/rms6_geostrophic_boundary_w0_day0560h00_two_layer_alternative.mat'),'omega','omega_QG','rhs','err_rms','lat');

end
toc