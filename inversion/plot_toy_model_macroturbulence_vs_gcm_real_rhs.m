% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (r-simulations) iteratively using Dirichlet-BC w=0 
% on the boundaries.

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
r_gcm = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

[lambda_vs_r,lambda_vs_r_aggregate,lambda_vs_r_max,lambda_vs_r_mean] = deal(zeros(10,1));

tic
for tt = 1:10

tt    
    
% Define Reduction Factor

r=r_gcm(tt); 

% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% inversion counter

lambda_vs_r_new = 0;
counter = 0;

% deformation radius

%load('figures_paper/deformation_radius_r_simulation.mat');
load('figures_paper/deformation_radius_r_simulation_p700_p400_average.mat');
Ld = Ld_vs_r(11,tt);

for i= 9:9 %9:12
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


% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta1 = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

err_rms = zeros(length(100),1);

lambda_vs_r_aggregate_new = 0;
lambda_vs_r_max_new = 0;
lambda_vs_r_mean_new = 0;
lambda_vs_lat = zeros(length(lat),1);
counter_max = 0;

w_aggregate = [];
for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);
sigma_accu(:,:,:,t) = -Ra * T(:,:,:,t)./ (Level .* theta) .* d_dp(theta, level);

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

% start the toy model inversions

[holder,p500] = min(abs(level-500*1e2));

lat0 = 45.5; % gcm domain doesn't start/end exactly at 25-65 lat;
S = mean(mean(mean(sigma_accu(:,:,p500,t),1),2),4);
S = squeeze(S);

% = [];


for jj = 1:length(lat)
    
% f0 = 2*Omega*sind(lat(jj));
% dp = 800*1e2;
% 
% Ld = sqrt(S)*dp/f0;
% Ld = Ld/(2*sqrt(2));

%Ld = 5.79*1e5;

% nondimensional version

L = 2*pi*R*cosd(lat(jj));
L = L/Ld;

N = length(lon);
x = linspace(0,L,N+1);
dx = L/N;    

rhs_500 = squeeze(rhs(jj,:,p500,t));
rhs_500 = rhs_500';
rhs_500 = rhs_500/max(abs(rhs_500));
rhs_500 = rhs_500-mean(rhs_500);
rhs_500 = -rhs_500; % go from pressure to z coordinates

[rhs_test] = generate_spectrum(abs(fft(rhs_500')).^2,L);

[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r,rhs_500,dx);

w_aggregate = cat(1,w_aggregate,w_final);

% % dimensional version
% 
% L = 2*pi*R*cosd(lat(jj));
% N = length(lon);
% x = linspace(0,L,N+1);
% dx = L/N;    
% rhs_500 = squeeze(rhs(jj,:,p500,t));
% rhs_500 = rhs_500';
% rhs_500 = -rhs_500; % go from pressure to z coordinates   
% 
% [w_final,asymmetry_parameter] = toy_model_dimensional_generated_spectrum(r,rhs_500,S,dp/2,f0,dx);

lambda_vs_lat(jj) = asymmetry_parameter;
lambda_vs_r_new = lambda_vs_r_new + asymmetry_parameter;

counter = counter + 1;


end

% lambda of the aggregated w

lambda_vs_r_aggregate_new = lambda_vs_r_aggregate_new + asymmetry(w_aggregate);

ind = find(islocalmax(lambda_vs_lat)==1);
lambda_vs_r_mean_new = lambda_vs_r_mean_new + mean(lambda_vs_lat(ind));

%lambda_vs_lat = sort(lambda_vs_lat,'descend');

lambda_vs_r_max_new = lambda_vs_r_max_new + max(lambda_vs_lat);


counter_max = counter_max + 1;
end

end

lambda_vs_r(tt) = lambda_vs_r_new/counter;
lambda_vs_r_aggregate(tt) = lambda_vs_r_aggregate_new/counter_max;
lambda_vs_r_max(tt) = lambda_vs_r_max_new/counter_max;
lambda_vs_r_mean(tt) = lambda_vs_r_mean_new/counter_max;

end

% calculate lambda from theory using rhs spectrum

lambda_vs_r_spectrum = zeros(length(r_gcm),1);
skew_rhs = zeros(length(r_gcm),1);

load('/disk7/mkohl/inversion/spec_rhs_gcm_r001.mat');
%load('figures_paper/deformation_radius_r_simulation_p700_p400_average.mat');

%Ld = 1e6/(2*sqrt(2));
%Ld = 5.79*1e5; % using S500 from model; see file plot_w_spec_r_simulation
%Ld = 5.21*1e5; % using weighted-sine method
%k_gcm = k*Ld;

draws = 100;

for ii = 1:length(r_gcm)
    
ii
    
[pholder,p500] = min(abs(level-50000));
%Ld = Ld_vs_r(p500,ii);
%Ld = 3.5*1e5;

Ld = Ld_vs_r(11,ii);


k_gcm = k*Ld;
    
spec_rhs_gcm = spec_500_r(ii,:); % select the r=0.01 run
spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);

R = 6371000.0;
%lat = 41.2336;
lat = 45.5; % gcm domain doesn't start/end exactly at 25-65 lat;
L = 2*pi*R*cosd(lat);
L = L/Ld;

N = length(spec_rhs_gcm);
x = linspace(0,L,N+1);
dx = L/N;

for jj = 1:draws
        
[rhs] = generate_spectrum(spec_rhs_gcm,L);
skew_rhs(ii) = skew_rhs(ii) + skewness(real(rhs));
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_gcm(ii),real(rhs),dx);


%[~,asymmetry_parameter] = toy_model(r_gcm(ii),3.3);
%
lambda_vs_r_spectrum(ii) = lambda_vs_r_spectrum(ii) + asymmetry_parameter;



end

skew_rhs(ii) = skew_rhs(ii)/draws;
lambda_vs_r_spectrum(ii) = lambda_vs_r_spectrum(ii)/draws;



% sample w-profile

% figure
% plot(x,w_final)
% xlabel('x')
% ylabel('w')
% xlim([x(1) x(end)])
% title(['r=0.01, $\lambda=$',num2str(round(asymmetry_parameter,2))])

end

clear('w_final');

% plot the asymmetries

% load GCM fields

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% % change interpreter to latex
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

figure
semilogx(r_gcm,lambda_turb_omega,'b'); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'b--'); hold on;
semilogx(r_gcm,lambda_vs_r,'k'); hold on;
semilogx(r_gcm,lambda_vs_r_aggregate,'k-.'); hold on;
semilogx(r_gcm,lambda_vs_r_mean,'k--'); hold on;
semilogx(r_gcm,lambda_vs_r_max,'k:'); hold on;
semilogx(r_gcm,lambda_vs_r_spectrum,'r')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('GCM','3-D Inversion','Theory True RHS (Aggregate)','Theory True RHS (Mean)','Theory True RHS (Peaks)','Theory True RHS (Max)','Theory RHS Spectrum','Location','SouthWest'); legend boxoff;
set(gca,'FontSize',12)

