% Calculates Lambda(r) for a Random RHS. 

clear;
close all;
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;
r =logspace(-2,0,30);
%r = 0.1;
runs = 30;
[asymmetry] = deal(zeros(length(runs), length(r)));

list = {'rescale0.01','rescale0.02','rescale0.05','rescale0.1','rescale0.2','rescale0.3',...                     
'rescale0.4','rescale0.6','rescale0.8','rescale1.0'};

strat_level = [19;14;13;14;14;13;12;18;12;18];

for tt = 1:length(r)

i=10;

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.1/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));

Omega = 7.2921e-5; % Omega = 7.2921e-5
R = 6371000.0;
Ra = 287.04; %287.04
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
event_timespan = [1:100];
level_indices = [5:30]; %5:30
level = plevels(level_indices);
p_up = 20000;
p_w = 5000;
p_up_index = find(level==p_up);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

T_in= ncread(path_temp,'temp_plevel', start, count);


[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));

for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        T(:, :, k, t) = T_in(:, :, k, t)';
    end
end
    
    
f0 = 2 * Omega * sin(lat_series_original(97) / 180 * 3.1415926);


[A, B, C, J, J1, J2, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
   
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs,sigma_r,del2_J, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level)));
[omega_QG] = deal(nan(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));


for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));


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

end

sigma = mean(sigma,2);

Reduction = reduction_factor(r(tt),level); 

% Inversion

for trials=1:runs
    
trials
        
rhs1 = randn(length(event_latspan), length(event_lonspan), length(level));

w = randn(length(lat),length(lon),length(level)); 

for s = 1:30
  
 for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level)
          if w(l,m,n)>=0
           sigma_r(l, m, n) = 0*sigma(n); 
          else 
           sigma_r(l, m, n) = (1-Reduction(n))*sigma(n);  
          end 
      end
    end
 end

J = w.*sigma_r;


for k=1:size(level)
    
    del2_J(:,:,k) = spherical_laplacian(J(:,:,k),phi,dphi,dlambda);
    
end
 
rhs_mod = rhs1+del2_J;

w = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                 sigma, dphi, dlambda, rhs_mod(:,:,:), omega_b(:,:,:));

end

asymmetry(trials) = Lambda(w(:,:,1:strat_level(tt)));

end

asymmetry_r(tt) = mean(asymmetry);

end