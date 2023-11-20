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

list = {'abs4.0_T170'};

strat_level_turbulent = [8,9,10,10,11,11,12, 12, 13, 13, 14, 14, 15,16,18];

tic
for tt = 1:1
tt
for i= 22:59
    i
err_rms = zeros(20,1);
% Define Path 

% Define Path 

path_temp = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/temp_day',num2str(5 *(521+i)),'h00.nc');
path_u = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/u_day',num2str(5 *(521+i)),'h00.nc');
path_v = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/v_day',num2str(5 *(521+i)),'h00.nc');
path_omega = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/omega_day',num2str(5 *(521+i)),'h00.nc');

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

event_latspan = [186:215]; %83 -111
event_lonspan = [1:512]; % 1:256
event_timespan = [1:20];
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
u_in= ncread(path_u,'u_plevel', start, count);
v_in= ncread(path_v,'v_plevel', start, count);
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

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

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

count = 0;

w_old = zeros(length(lat),length(lon),length(level_place));

% w_new = randn(size(w_old));

% w_new = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%   sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
sigma_place(:), dphi, dlambda, rhs_place(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));

diff = w_new-w_old;

% Pre-Allocation

sigma_r = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));
del2_J = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));

%while rms(diff(:))>=10^(-4) && count<400

%plot(squeeze(w_new(15,:,10))); hold on;
    
while rms(diff(:))>=1*10^(-4) && count<50
    
count = count +1;
       
 for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level_place)
          if w_new(l,m,n)>=0
           sigma_r(l, m, n) = 0*sigma_place(n); 
          else 
           sigma_r(l, m, n) = (1-r_place(l,m,n))*sigma_place(n);  
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

% w_new = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%                       sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));
                   
w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
 sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));               

diff = w_new-w_old;

% format long
% error(count) = rms(diff(:))
% asym(count) = Lambda_weighted(w_new(11:26,:,1:strat_level_turbulent(tt),1),lat(11:26));
% 
% plot(squeeze(w_new(15,:,10))); hold on;
                                                                                                                  
end
err_rms(t) = rms(diff(:));
omega_QG(:,:,bound:end,t) = w_new;
end

% figure(6)
% plot(lon,omega(11,:,7+bound,1)); hold on;
% plot(lon,w_new(11,:,7,1)); hold on;
% plot(lon,test(11,:,7,1)); hold on;
% legend('omega','QG','dry')
% xlabel('lon')
% ylabel('w(Pa/s)')
% title('w at 500hPA and 40 lat, r=0.4')

save(strcat('/net/graupel/archive1/users/mkohl/inversion/output/abs_simulation/',list{tt},'/rms1em4_geostrophic_boundary_w_day',num2str(5 *(521+i)),'_T170.mat'),'omega','omega_QG','r','rhs','lat','err_rms');
end

end
toc