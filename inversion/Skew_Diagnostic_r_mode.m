% Diagnose the Skewness of different terms on the RHS of omega-equation.

clear;

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
    
    
f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, J, Qx, Qy, zeta, ...
    du_dlambda, dv_dlambda,du_dp,dv_dp, du_dphi, dv_dphi, diab ...
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
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);

for k = 1:length(level)
    
    zeta(:,:,k,t) = curl(phi,u(:,:,k,t),v(:,:,k,t),dphi,dlambda);

end

du_dp(:,:,:,t) = d_dp(u(:,:,:,t),level);
dv_dp(:,:,:,t) = d_dp(v(:,:,:,t),level);

%rr = r_parameter1(T(:,:,:,t),Level,level);
rr = ones(size(T(:,:,:,t)));
rr(omega(:,:,:,t)<=0) = r;
diab(:,:,:,t) = (1-rr).*omega(:,:,:,t).*S(:,:,:,t);
diab(:,:,:,t) = diab(:,:,:,t).*(1e5./Level).^kappa;

for k = 1:length(level)
diab(:,:,k,t) =  spherical_laplacian(diab(:,:,k,t), phi, dphi, dlambda);   
end

end

up = u - mean(u,2);
vp = v - mean(v,2);
zetap = zeta-mean(zeta,2);

um = repmat(mean(u,2),1,size(u,2),1,1);
vm = repmat(mean(v,2),1,size(u,2),1,1);
zetam = repmat(mean(zeta,2),1,size(u,2),1,1);

[A, A1, A2, A3] = traditional_A(phi,level,u,v,zeta,event_timespan,dphi,dlambda,f0,beta);

[Amp, A1mp, A2mp, A3mp] = traditional_A(phi,level,um,vm,zetap,event_timespan,dphi,dlambda,f0,beta);

[App, A1pp, A2pp, A3pp] = traditional_A(phi,level,up,vp,zetap,event_timespan,dphi,dlambda,f0,beta);


B = traditional_B(phi,level,u,v,T,event_timespan,dphi,dlambda);


Deform = Deformation(phi,level, u, v, T, event_timespan, dphi, dlambda);

rhs = 2*A1+A3+Deform;

%Deformation = A+B - 2*A1-A3;

Product = (du_dp+dv_dp).*zeta;

skew_rhs(tt) = -Skewness(rhs(:,:,index(tt),:));

skew_diab(tt) = -Skewness(diab(:,:,index(tt),:));

skew_A1(tt) = -Skewness(A1(:,:,index(tt),:));

skew_A1mp(tt) = -Skewness(A1mp(:,:,index(tt),:));

skew_A1pp(tt) = -Skewness(A1pp(:,:,index(tt),:));

skew_A3(tt) = -Skewness(A3(:,:,index(tt),:));

skew_Deform(tt) = -Skewness(Deform(:,:,index(tt),:));

skew_shear_u(tt) = -Skewness(du_dp(:,:,index(tt),:));

skew_shear_v(tt) = -Skewness(dv_dp(:,:,index(tt),:));

skew_vor(tt) = -Skewness(zeta(:,:,index(tt),:));

skew_Product(tt) = -Skewness(Product(:,:,index(tt),:));

% Norm

rhs = squeeze(rhs(:,:,index(tt),:)); 
diab = squeeze(diab(:,:,index(tt),:));
A1 = squeeze(A1(:,:,index(tt),:));
A1mp = squeeze(A1mp(:,:,index(tt),:));
A1pp = squeeze(A1pp(:,:,index(tt),:));
A3 = squeeze(A3(:,:,index(tt),:));
Deform = squeeze(Deform(:,:,index(tt),:));

norm_rhs(tt) = norm(rhs(:));
norm_diab(tt) = norm(diab(:));
norm_A1(tt) = norm(A1(:));
norm_A1mp(tt) = norm(A1mp(:));
norm_A1pp(tt) = norm(A1pp(:));
norm_A3(tt) = norm(A3(:));
norm_Deform(tt) = norm(Deform(:));

end

% Make a plot of the w-profile in ascending regions

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

w_up_pp  = d_d2p(w_up,level);

coeff = (500*1e2)^2/2*w_up_pp./w_up;

figure(100)
plot(w_up(5:end),level(5:end)/100)
xlabel('w_{up}')
ylabel('p(hPa)')
title('r=0.01 mode')
set(gca,'Ydir','reverse');


figure(1)
semilogx(R_factor,skew_rhs); hold on;
semilogx(R_factor,skew_diab); hold on;
semilogx(R_factor,skew_A1); hold on;
semilogx(R_factor,skew_A1mp); hold on;
semilogx(R_factor,skew_A1pp); hold on;
semilogx(R_factor,skew_A3); hold on;
semilogx(R_factor,skew_Deform)
xlabel('r')
ylabel('Skew')
title('Skew at 500hPA for r-mode')
%legend('RHS','diab','vor-adv','beta','Deformation','location','northeast')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
%savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_r_mode_1.fig')

figure(2)
loglog(R_factor,norm_rhs); hold on;
loglog(R_factor,norm_diab); hold on;
loglog(R_factor,norm_A1); hold on;
loglog(R_factor,norm_A1mp); hold on;
loglog(R_factor,norm_A1pp); hold on;
loglog(R_factor,norm_A3); hold on;
loglog(R_factor,norm_Deform)
xlabel('r')
ylabel('Norm')
title('Norm at 500hPA for r-mode')
%legend('RHS','diab','Vor Adv','Beta','Deformation','location','northeast')
%legend('RHS','diab','Vor Adv','Beta','Deformation','location','northeast')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)


% figure(2)
% %set(gcf,'DefaultAxesFontSize', 14)
% semilogx(R_factor,Skew_total); hold on;
% semilogx(R_factor,Skew_shear_u); hold on;
% semilogx(R_factor,Skew_shear_v); hold on;
% semilogx(R_factor,Skew_vor); hold on;
% semilogx(R_factor,Skew_Product)
% xlabel('r')
% ylabel('Skew')
% title('Skew J(\tau,del \phi) Decomposition at 500hPA for r-mode')
% legend('Total','du/dp','dv/dp','zeta','du/dp*zeta','location','northeast')
% legend boxoff
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
%savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_r_mode_2.fig')

% figure(1)
% semilogx(R_factor,Skew);
% xlabel('r')
% ylabel('Skew')
% title('Skew J(tau,del tau) at 500hPA for Mode r-simulation')
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_r_mode.fig')
% 
% figure(2)
% semilogx(R_factor,Lambda_parameter)
% xlabel('r')
% ylabel('\lambda')
% title('Lambda J(tau,del tau) at 500hPA for Mode r-simulation')
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/lambda_r_mode.fig')
% 
% 



