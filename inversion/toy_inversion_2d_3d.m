% Compare 2d vs 3d inversions

clear; close all;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

%strat_level = [19;14;13;14;14;13;12;18;12;18];

strat_level = 14*ones(size(R_factor));

index = [16,11,11,11,11,11,11,16,11,16];

index_up = index+4;
index_low = index-4;

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

f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);

    
%f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
%f = f0*ones(size(lat));
%beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
% = beta*ones(size(lat));  

[A, B, C, J, Qx, Qy, zeta, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, diab...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed,omega_adv,omega_diab] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
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
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);

for k = 1:length(level)
    
    zeta(:,:,k,t) = curl(phi,u(:,:,k,t),v(:,:,k,t),dphi,dlambda);

end

r = reduction_factor(R_factor(tt),level);

r = repmat(reshape(r, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);

%[omega_adv(:,:,:,t), omega_diab(:,:,:,t)] = inversion(u(:,:,:,t),v(:,:,:,t),omega(:,:,:,t),T(:,:,:,t),r,sigma(:,t),f,beta,level,Phi,phi,dphi,lambda,dlambda);

[omega_adv(:,:,:,t), omega_diab(:,:,:,t)] = inversion_toy_rhs(u(:,:,:,t),v(:,:,:,t),omega(:,:,:,t),T(:,:,:,t),r,sigma(:,t),f,beta,level,Phi,phi,dphi,lambda,dlambda);

end

end


% rhs = randn(size(omega));

% w_new = SIP_diagnostic(R, f0, size(omega,2), size(omega,1), length(level), phi, lambda, level, ...
% sigma(:), dphi, dlambda, rhs(:,:,:),  omega_b(:,:,:));



%[omega_adv, omega_diab] = inversion_toy_rhs(u,v,omega,T,r,sigma,f,beta,level,Phi,phi,dphi,lambda,dlambda);

%[omega_adv, omega_diab] = inversion_constant_f(u,v,omega,T,r,sigma,f0,beta,level,Phi,phi,dphi,lambda,dlambda);

omega = omega-mean(omega,2);
omega_diab = omega_diab - mean(omega_diab,2);


skew_omega_gcm(tt) = mean(-Skewness(omega(:,:,1:17,:)));

lambda_omega_gcm(tt) = mean(Lambda(omega(:,:,1:17,:)));

skew_omega_qg(tt) = mean(-Skewness(omega_diab(:,:,1:17,:)));

lambda_omega_qg(tt) = mean(Lambda(omega_diab(:,:,1:17,:)));

skew_omega_gcm_500(tt) = mean(-Skewness(omega(:,:,11,:)));

lambda_omega_gcm_500(tt) = mean(Lambda(omega(:,:,11,:)));

skew_omega_qg_500(tt) = mean(-Skewness(omega_diab(:,:,11,:)));

lambda_omega_qg_500(tt) = mean(Lambda(omega_diab(:,:,11,:)));


end
toc

figure(1)
semilogx(R_factor,lambda_omega_gcm); hold on;
semilogx(R_factor,lambda_omega_qg);
xlabel('r')
ylabel('\lambda')
legend('GCM','QG'); legend boxoff
saveas(gcf,'/net/halo/disk28/disk7/mkohl/inversion/rhs_toy_model/lambda_gcm_vs_qg_rhs500_const','epsc')

figure(2)
semilogx(R_factor,lambda_omega_gcm_500); hold on;
semilogx(R_factor,lambda_omega_qg_500);
xlabel('r')
ylabel('\lambda 500')
legend('GCM','QG'); legend boxoff
saveas(gcf,'/net/halo/disk28/disk7/mkohl/inversion/rhs_toy_model/lambda_500_gcm_vs_qg_rhs500_const','epsc')

