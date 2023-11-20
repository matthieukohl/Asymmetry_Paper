% plot spectrum of w

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

[wavenumber_up_vs_r,wavenumber_down_vs_r,wavenumber_total_vs_r] = deal(zeros(10,1));

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

wavenumber_total = 0;
wavenumber_up = 0;
wavenumber_down = 0;
counter = 0;

for i= 9:12 %9:12
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/ug_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/vg_day0',num2str(25 * i, '%.3d'),'h00.nc');
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
%event_latspan = [94:98]; %
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
    
% % varying f    
% f = 2 * Omega * sind(lat);
% beta = 1 / R * d_dphi(f, dphi);

% constant f0
lat0 = lat(1)+(lat(end)-lat(1))/2;
f = 2 * Omega * sind(lat0*ones(size(lat)));
beta = 2 * Omega * cosd(lat0*ones(size(lat)))/R;

% calculate the rhs

[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

err_rms = zeros(length(100),1);
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

% calculate the laplacian

[a,p500] = min(abs(level-500*1e2));

% now evaluate mean(lap(w)^2)/mean(w^2)
% over the whole domain/ over ascending regions
% over descending regions

for t = 1:length(event_timespan)
    for kk = p500:p500

rhs_2d = rhs(:,:,kk,t);
%w_2d = w_2d - mean(mean(w_2d));
lap_rhs = spherical_laplacian(rhs_2d, phi, dphi, dlambda);


% paul's formulation

% square formulation

% wavenumber_total = wavenumber_total + mean(mean(lap_w.^2))/(mean(mean(w_2d.^2)));
% wavenumber_up = wavenumber_up + mean(lap_w(w_2d<0).^2)/(mean(w_2d(w_2d<0).^2));
% wavenumber_down = wavenumber_down + mean(lap_w(w_2d>0).^2)/(mean(w_2d(w_2d>0).^2));

% absolute value formulation

wavenumber_total = wavenumber_total + mean(mean(abs(lap_rhs)))/(mean(mean(abs(rhs_2d))));
wavenumber_up = wavenumber_up + mean(mean(lap_rhs(rhs_2d>0)))/(mean(mean(rhs_2d(rhs_2d>0))));
wavenumber_down = wavenumber_down + mean(mean(lap_rhs(rhs_2d<0)))/(mean(mean(rhs_2d(rhs_2d<0))));

% wavenumber_total = wavenumber_total + mean(mean(abs(lap_w)))/(mean(mean(abs(w_2d))));
% wavenumber_up = wavenumber_up + mean(mean(abs(lap_w(w_2d<0))))/(mean(mean(abs(w_2d(w_2d<0)))));
% wavenumber_down = wavenumber_down + mean(mean(abs(lap_w(w_2d>0))))/(mean(mean(abs(w_2d(w_2d>0)))));



% my formulation

% % absolute value

%ratio = lap_w./w_2d;

% wavenumber_total = wavenumber_total + mean(mean(abs(ratio)));
% wavenumber_up = wavenumber_up + mean(mean(abs(ratio(w_2d<0))));
% wavenumber_down = wavenumber_down + mean(mean(abs(ratio(w_2d>0))));

% wavenumber_total = wavenumber_total + mean(mean(ratio.^2));
% wavenumber_up = wavenumber_up + mean(mean(ratio(w_2d<0).^2));
% wavenumber_down = wavenumber_down + mean(mean(ratio(w_2d>0).^2));


counter = counter + 1;

    end
end


end

wavenumber_total = wavenumber_total/counter;
wavenumber_up = wavenumber_up/counter;
wavenumber_down = wavenumber_down/counter;

% wavenumber_total_vs_r(tt) = (wavenumber_total/4)^(1/4);
% wavenumber_up_vs_r(tt) = (wavenumber_up/4)^(1/4);
% wavenumber_down_vs_r(tt) = (wavenumber_down/4)^(1/4);

wavenumber_total_vs_r(tt) = (wavenumber_total/2)^(1/2);
wavenumber_up_vs_r(tt) = (-wavenumber_up/2)^(1/2);
wavenumber_down_vs_r(tt) = (-wavenumber_down/2)^(1/2);

% wavenumber_total_vs_r(tt) = (wavenumber_total/2)^(1/2);
% wavenumber_up_vs_r(tt) = (wavenumber_up/2)^(1/2);
% wavenumber_down_vs_r(tt) = (wavenumber_down/2)^(1/2);

% wavenumber_total_vs_r(tt) = (wavenumber_total/2)^(-1/2);
% wavenumber_up_vs_r(tt) = (wavenumber_up/2)^(-1/2);
% wavenumber_down_vs_r(tt) = (wavenumber_down/2)^(-1/2);


end

Ld = 5.79*1e5;
wavenumber_up_vs_r = wavenumber_up_vs_r*Ld;
wavenumber_down_vs_r = wavenumber_down_vs_r*Ld;
wavenumber_total_vs_r = wavenumber_total_vs_r*Ld;

% estimate total k from 2/k = 1/ku + 1/kd relation
% k = 2*ku*kd/(ku+kd)

wavenumber_total_estimate_vs_r = 2*wavenumber_up_vs_r.*wavenumber_down_vs_r./(wavenumber_up_vs_r+wavenumber_down_vs_r);

set(groot,'DefaultTextInterpreter','default')
set(groot,'DefaultLegendInterpreter','default')

figure
semilogx(R_factor,wavenumber_down_vs_r,'b','Linewidth',1.4); hold on;
semilogx(R_factor,wavenumber_up_vs_r,'r','Linewidth',1.4); hold on;
semilogx(R_factor,wavenumber_total_vs_r,'k','Linewidth',1.4); hold on;
semilogx(R_factor,wavenumber_total_estimate_vs_r,'k--','Linewidth',1.4);
xlabel('Reduction factor r')
ylabel('Wavenumber k')
%title('\rm mean((Lap w)^2) / mean(w^2)')
title('\rm square formulation')
set(gca,'FontSize',12);
legend('k_d','k_u','k','k estimate'); legend boxoff; 