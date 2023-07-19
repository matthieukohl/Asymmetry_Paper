% Diagnose the Skewness of different terms on the RHS of omega-equation.

clear; close all

% load the macroturbulent runs

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

R_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];


index = [16,11,11,11,11,11,11,16,11,16];

index_up = index+4;
index_low = index-4;

for tt = 1:1

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

for i= 9:9 %9:12
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/u_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/v_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');

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
%event_latspan = [93:101]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [1:100];
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
    
f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, J, Qx, Qy, zeta, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, diab...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);


end

end

end

%omega = omega-mean(omega,2);

w_up  = zeros(size(omega));

for ii = 1:size(omega,1)
    for jj = 1:size(omega,2)
        for kk = 1:size(omega,3)
            for tt = 1:size(omega,4)
                
 if omega(ii,jj,kk,tt)<0
     w_up(ii,jj,kk,tt) = omega(ii,jj,kk,tt);
 else
     w_up(ii,jj,kk,tt) = NaN;
 end
            end
        end
    end
end

w_up = nanmean(nanmean(nanmean(w_up,1),2),4);
w_up = squeeze(w_up);

level_macro = level;
w_up_macro = w_up/(max(abs(w_up)));

% load the mode runs

list = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};
 
R_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

index = [11,11,11,11,11,11,11,11,11,11];

for tt = 1:1
    
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
    
end


w_up  = zeros(size(omega));

%omega = omega-mean(omega,2);

for ii = 1:size(omega,1)
    for jj = 1:size(omega,2)
        for kk = 1:size(omega,3)
            for tt = 1:size(omega,4)
                
 if omega(ii,jj,kk,tt)<0
     w_up(ii,jj,kk,tt) = omega(ii,jj,kk,tt);
 else
     w_up(ii,jj,kk,tt) = NaN;
 end
            end
        end
    end
end

w_up = nanmean(nanmean(nanmean(w_up,1),2),4);
w_up = squeeze(w_up);

level_mode = level;
w_up_mode = w_up/(max(abs(w_up)));

% make a plot
r = 0.01;
f = r + 0.5*(1-r)*(1-tanh((level-20000)/5000));

figure(1)
plot(-w_up_mode,level_mode/100); hold on;
plot(-w_up_macro,level_macro/100); hold on
plot(-f,level/100);
ylabel('p(hPa)')
xlabel('norm. magnitude')
legend('mode','macro','-r(p)','Location','north'); legend boxoff
set(gca,'Ydir','reverse')
title('w-profile')


S = nanmean(nanmean(nanmean(sigma_accu,1),2),4);
S = squeeze(S);
diab_tend_mode = d_dp((1-f).*S.*w_up_mode,level);

diab_tend_macro = d_dp((1-f).*S.*w_up_macro,level);



figure(2)
plot(diab_tend_mode,level_mode/100); hold on;
plot(diab_tend_macro,level_macro/100); hold on
ylabel('p(hPa)')
xlabel('magnitude')
title('d_{dp}(1-r)S\omega_{norm}')
legend('mode','macro','Location','northwest'); legend boxoff
set(gca,'Ydir','reverse')

