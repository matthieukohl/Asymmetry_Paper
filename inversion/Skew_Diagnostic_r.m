% Diagnose the Skewness of different terms on the RHS of omega-equation.

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

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

for i= 12:12 %9:12
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
S(:,:,:,t) = -T(:,:,:,t)./(theta).*d_dp(theta,level);

for k = 1:length(level)
    
    zeta(:,:,k,t) = curl(phi,u(:,:,k,t),v(:,:,k,t),dphi,dlambda);

end

%rr = r_parameter1(T(:,:,:,t),Level,level);
rr = ones(size(T(:,:,:,t)));
rr(omega(:,:,:,t)<=0) = r;
diab(:,:,:,t) = (1-rr).*omega(:,:,:,t).*S(:,:,:,t);
diab(:,:,:,t) = diab(:,:,:,t).*(1e5./Level).^kappa;

for k = 1:length(level)
diab(:,:,k,t) =  spherical_laplacian(diab(:,:,k,t), phi, dphi, dlambda);   
end

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

skew_omega(tt) = -Skewness(omega(:,:,index(tt),:));

lambda_omega(tt) = Lambda(omega(:,:,index(tt),:));

skew_rhs(tt) = -Skewness(rhs(:,:,index(tt),:));

skew_diab(tt) = -Skewness(diab(:,:,index(tt),:));

skew_A1(tt) = -Skewness(A1(:,:,index(tt),:));

skew_A1mp(tt) = -Skewness(A1mp(:,:,index(tt),:));

skew_A1pp(tt) = -Skewness(A1pp(:,:,index(tt),:));

skew_A3(tt) = -Skewness(A3(:,:,index(tt),:));

skew_A3mp(tt) = -Skewness(A3mp(:,:,index(tt),:));

skew_A3pp(tt) = -Skewness(A3pp(:,:,index(tt),:));

skew_Deform(tt) = -Skewness(Deform(:,:,index(tt),:));


% Norm

u500 = squeeze(u(:,:,index(tt),:));
v500 = squeeze(v(:,:,index(tt),:));
u300 = squeeze(u(:,:,index_up(tt),:));
v300 = squeeze(v(:,:,index_up(tt),:));
u700 = squeeze(u(:,:,index_low(tt),:));
v700 = squeeze(v(:,:,index_low(tt),:));

w500 = squeeze(omega(:,:,index(tt),:));
rhs = squeeze(rhs(:,:,index(tt),:)); 
diab = squeeze(diab(:,:,index(tt),:));
A1 = squeeze(A1(:,:,index(tt),:));
A1mp = squeeze(A1mp(:,:,index(tt),:));
A1pp = squeeze(A1pp(:,:,index(tt),:));
A3 = squeeze(A3(:,:,index(tt),:));
A3mp = squeeze(A3mp(:,:,index(tt),:));
A3pp = squeeze(A3pp(:,:,index(tt),:));
Deform = squeeze(Deform(:,:,index(tt),:));

norm_rhs(tt) = norm(rhs(:));
norm_diab(tt) = norm(diab(:));
norm_A1(tt) = norm(A1(:));
norm_A1mp(tt) = norm(A1mp(:));
norm_A1pp(tt) = norm(A1pp(:));
norm_A3(tt) = norm(A3(:));
norm_A3mp(tt) = norm(A3mp(:));
norm_A3pp(tt) = norm(A3pp(:));
norm_Deform(tt) = norm(Deform(:));

RHS = rhs(:,:,end);
mean_rhs(tt) = mean(RHS(:));
std_rhs(tt) = std(RHS(:));
skew_rhs(tt) = skewness(RHS(:));
kurt_rhs(tt) = kurtosis(RHS(:));


% u500 = u500-mean(u500,2);
% v500 = v500-mean(v500,2);
% u700 = u700-mean(u700,2);
% v700 = v700-mean(v700,2);
% u300 = u300-mean(u300,2);
% v300 = v300-mean(v300,2);
% w500 = w500-mean(w500,2);
% u_phi = 1/2*(u300+u700);
% u_tau = 1/2*(u300-u700);
% rhs = rhs-mean(rhs,2);
% En_uv = 1/2*(u500.^2 + v500.^2);
% 
% En_phi = 1/2*((0.5*u300+0.5*u700).^2+(0.5*v300+0.5*v700).^2);
% 
% En_phi = En_phi - mean(En_phi,2);
% En_phi = En_phi - mean(En_phi,1);
% 
% En_tau = 1/2*((0.5*u300-0.5*u700).^2+(0.5*v300-0.5*v700).^2);
% 
% En_tau = En_tau - mean(En_tau,2);
% En_tau = En_tau - mean(En_tau,1);
% 
% En_w = 1/2*(w500.^2);
% % En_w = En_w - mean(En_w,1);
% 
% 
% power_rhs = 0;
% 
% counter = 0;
% for ii = 1:size(rhs,1)
%     for t = 1:size(rhs,3)
% rhsf = fft(rhs(ii,:,t));
% power = abs(rhsf).^2;
% power_rhs = power_rhs + power;
% counter = counter +1;
%     end
% end
% power_rhs = power_rhs/(size(rhs,1)*size(rhs,3));
% 
% 
% power_rhs_ky = 0;
% 
% 
% counter = 0;
% for ii = 1:size(rhs,1)
%     for t = 1:size(rhs,3)
% rhsf = fft(rhs(ii,:,t));
% power = abs(rhsf).^2;
% power_rhs_ky = power_rhs_ky + power;
% counter = counter +1;
%     end
% end
% power_rhs_ky = power_rhs_ky/(size(rhs,2)*size(rhs,3));
% 
% [minValue,latp] = min(abs(lat-41));
% 
% 
% ny = size(rhs,1);
% Ly = 2*pi*R*(abs(lat(end)-lat(1)))/360;
% ky = 2*pi/Ly*[0:ny,-ny/2+1:-1];
% ky = ky*1e6/sqrt(2);
% 
% nx = size(rhs,2);
% L = 2*pi*R*cosd(lat(latp));
% k = 2*pi/L*[0:nx,-nx/2+1:-1];
% k = k*1e6/sqrt(2);
% n = [0:nx,-nx/2+1:-1];
% 
% 
% kx_centroid(tt) = sum(k(1:nx/2).*power_rhs(1:nx/2))./(sum(power_rhs(1:nx/2)));
% 
% ky_centroid(tt) = sum(ky(1:round(ny/2)).*power_rhs_ky(1:round(ny/2)))./(sum(power_rhs_ky(1:round(ny/2))));
% 
% 
% power_w = 0;
% 
% counter = 0;
% for ii = 1:size(rhs,1)
%     for t = 1:size(rhs,3)
% rhsf = fft(w500(ii,:,t));
% power = abs(rhsf).^2;
% power_w = power_w + power;
% counter = counter +1;
%     end
% end
% power_w = power_w/(size(rhs,1)*size(rhs,3));
% 
% 
% power_w_ky = 0;
% 
% 
% counter = 0;
% for ii = 1:size(rhs,1)
%     for t = 1:size(rhs,3)
% rhsf = fft(w500(ii,:,t));
% power = abs(rhsf).^2;
% power_w_ky = power_w_ky + power;
% counter = counter +1;
%     end
% end
% power_w_ky = power_w_ky/(size(rhs,2)*size(rhs,3));
% 
% 
% kx_w_centroid(tt) = sum(k(1:nx/2).*power_w(1:nx/2))./(sum(power_w(1:nx/2)));
% 
% ky_w_centroid(tt) = sum(ky(1:round(ny/2)).*power_w_ky(1:round(ny/2)))./(sum(power_w_ky(1:round(ny/2))));


end

figure(200)
RHS = rhs(:,:,end);
% %[ minA , maxA ] = bounds(rhs(:));
% BE = linspace(-4e-16,4e16,1e5);
% 
% for tt = 1:size(rhs,3)
%   counts(:,tt) = histcounts(rhs(:,:,tt),BE);
% end
% counts_av = mean(counts,2);
% 
% histogram(counts_av,'BinEdges',BE);

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
title('r=0.01 macroturbulent')
set(gca,'Ydir','reverse');

figure(200)
plot(coeff(5:end),level(5:end)/100)
xlabel('w_{pp}/w')
ylabel('p(hPa)')
set(gca,'Ydir','reverse');


histogram(RHS,1000,'Normalization','pdf');
skew = skewness(RHS(:));
kurt = kurtosis(RHS(:));
stdd = std(RHS(:));
title(['st=',num2str(round(stdd,1)),'  skew=',num2str(round(skew,1)),'  kurt=',num2str(round(kurt,1))])




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
legend('RHS','diab','vor-adv','beta','Deformation','location','northeast')
legend boxoff
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% %savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_r_mode_1.fig')
% 
figure(2)
loglog(R_factor,norm_rhs); hold on;
loglog(R_factor,norm_diab); hold on;
loglog(R_factor,norm_A1); hold on;
loglog(R_factor,norm_A1mp); hold on
loglog(R_factor,norm_A1pp); hold on;
loglog(R_factor,norm_A3); hold on;
loglog(R_factor,norm_A3mp); hold on;
loglog(R_factor,norm_A3pp); hold on;
loglog(R_factor,norm_Deform)
xlabel('r')
ylabel('Norm')
title('Norm at 500hPA for r-mode')
legend('RHS','diab','Vor Adv','Beta','Deformation','location','northeast')
%legend('RHS','Vor Adv','Beta','Deformation','location','northeast')
legend boxoff
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)



% figure(1)
% semilogx(R_factor,Skew);
% xlabel('r')
% ylabel('Skew')
% title('Skew J(tau,del tau) at 500hPA for turbulent r-simulation')
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_r_turbulent.fig')
% 
% figure(2)
% semilogx(R_factor,Lambda_parameter)
% xlabel('r')
% ylabel('\lambda')
% title('Lambda J(tau,del tau) at 500hPA for turbulent r-simulation')
% savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/lambda_r_turbulent.fig')
% 
u500 = u500-mean(u500,2);
v500 = v500-mean(v500,2);
u700 = u700-mean(u700,2);
v700 = v700-mean(v700,2);
u300 = u300-mean(u300,2);
v300 = v300-mean(v300,2);
w500 = w500-mean(w500,2);
u_phi = 1/2*(u300+u700);
u_tau = 1/2*(u300-u700);
rhs = rhs-mean(rhs,2);
En_uv = 1/2*(u500.^2 + v500.^2);

En_phi = 1/2*((0.5*u300+0.5*u700).^2+(0.5*v300+0.5*v700).^2);

En_phi = En_phi - mean(En_phi,2);
En_phi = En_phi - mean(En_phi,1);

En_tau = 1/2*((0.5*u300-0.5*u700).^2+(0.5*v300-0.5*v700).^2);

En_tau = En_tau - mean(En_tau,2);
En_tau = En_tau - mean(En_tau,1);

En_w = 1/2*(w500.^2);


% plot spectrum of the rhs at 500hPa
[Lat,Lon] = meshgrid(lon,lat);
%latp = 12;
[minValue,latp] = min(abs(lat-41));


power_rhs = 0;

counter = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
rhsf = fft(rhs(ii,:,t));
power = abs(rhsf).^2;
power_rhs = power_rhs + power;
counter = counter +1;
    end
end
power_rhs = power_rhs/(size(rhs,1)*size(rhs,3));


power_rhs_ky = 0;


counter = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
rhsf = fft(rhs(ii,:,t));
power = abs(rhsf).^2;
power_rhs_ky = power_rhs_ky + power;
counter = counter +1;
    end
end
power_rhs_ky = power_rhs_ky/(size(rhs,2)*size(rhs,3));

power_u_tau_kx = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(u_tau(ii,:,end));
power = abs(ff).^2;
power_u_tau_kx = power_u_tau_kx + power;
    end
end
power_u_tau_kx = power_u_tau_kx/(size(rhs,1)*size(rhs,3));

power_u_tau_ky = 0;
for ii = 1:size(rhs,2)
    for t = 1:size(rhs,3)
ff = fft(u_tau(:,ii,end));
power = abs(ff).^2;
power_u_tau_ky = power_u_tau_ky + power;
    end
end
power_u_tau_ky = power_u_tau_ky/(size(rhs,2)*size(rhs,3));

power_u_phi_kx = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(u_phi(ii,:,end));
power = abs(ff).^2;
power_u_phi_kx = power_u_phi_kx + power;
    end
end
power_u_phi_kx = power_u_phi_kx/(size(rhs,1)*size(rhs,3));

power_u_phi_ky = 0;
for ii = 1:size(rhs,2)
    for t = 1:size(rhs,3)
ff = fft(u_phi(:,ii,end));
power = abs(ff).^2;
power_u_phi_ky = power_u_phi_ky + power;
    end
end
power_u_phi_ky = power_u_phi_ky/(size(rhs,2)*size(rhs,3));


power_uv = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(En_uv(ii,:,t));
power = abs(ff).^2;
power_uv = power_uv + power;
    end
end
power_uv = power_uv/(size(rhs,1)*size(rhs,3));

power_phi = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(En_phi(ii,:,t));
power = abs(ff).^2;
power_phi = power_phi + power;
    end
end
power_phi = power_phi/(size(rhs,1)*size(rhs,3));

power_phi_ky = 0;
for ii = 1:size(rhs,2)
    for t = 1:size(rhs,3)
ff = fft(En_phi(:,ii,t));
power = abs(ff).^2;
power_phi_ky = power_phi_ky + power;
    end
end
power_phi_ky = power_phi_ky/(size(rhs,2)*size(rhs,3));

power_tau = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(En_tau(ii,:,end));
power = abs(ff).^2;
power_tau = power_tau + power;
    end
end
power_tau = power_tau/(size(rhs,1)*size(rhs,3));

power_tau_ky = 0;

for ii = 1:size(rhs,2)
    for t = 1:size(rhs,3)
ff = fft(En_tau(:,ii,end));
power = abs(ff).^2;
power_tau_ky = power_tau_ky + power;
    end
end
power_tau_ky = power_tau_ky/(size(rhs,2)*size(rhs,3));


power_w = 0;
for ii = 1:size(rhs,1)
    for t = 1:size(rhs,3)
ff = fft(En_w(ii,:,t));
power = abs(ff).^2;
power_w = power_w + power;
    end
end
power_w = power_w/(size(rhs,1)*size(rhs,3));

nx = size(rhs,2);
L = 2*pi*R*cosd(lat(latp));
k = 2*pi/L*[0:nx,-nx/2+1:-1];
k = k*1e6/sqrt(2);
n = [0:nx,-nx/2+1:-1];

figure(200)
RHS = rhs(:,:,end);
% %[ minA , maxA ] = bounds(rhs(:));
% BE = linspace(-4e-16,4e16,1e5);
% 
% for tt = 1:size(rhs,3)
%   counts(:,tt) = histcounts(rhs(:,:,tt),BE);
% end
% counts_av = mean(counts,2);
% 
% histogram(counts_av,'BinEdges',BE);


histogram(RHS,1000,'Normalization','pdf');
skew = skewness(RHS(:));
kurt = kurtosis(RHS(:));
stdd = std(RHS(:));
title(['st=',num2str(round(stdd,1)),'  skew=',num2str(round(skew,1)),'  kurt=',num2str(round(kurt,1))])

figure(1)
loglog(n(1:nx/2),power_rhs(1:nx/2))
ylabel('abs(fft(rhs))^2')
xlabel('n')
xlim([n(1) n(nx/2)])

figure('Position',[10 10 800 200])

subplot(1,2,1)

plot(k(1:nx/2),power_rhs(1:nx/2)/max(power_rhs(1:nx/2))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('kx')
xlim([k(1) k(nx/2)])
title('X Spectra RHS')


ny = size(rhs,1);
Ly = 2*pi*R*(abs(lat(end)-lat(1)))/360;
ky = 2*pi/Ly*[0:ny,-ny/2+1:-1];
ky = ky*1e6/sqrt(2);

subplot(1,2,2)

plot(ky(1:round(ny/2)),power_rhs(1:round(ny/2))/max(power_rhs(1:round(ny/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('ky')
xlim([ky(1) ky(round(ny/2))])
title('Y Spectra RHS')

figure('Position',[10 10 1200 500])


subplot(2,2,2)
plot(k(1:round(nx/2)),power_u_tau_kx(1:round(nx/2))/max(power_u_tau_kx(1:round(nx/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('kx')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra u_{\tau}')

subplot(2,2,4)
plot(ky(1:round(ny/2)),power_u_tau_ky(1:round(ny/2))/max(power_u_tau_ky(1:round(ny/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('ky')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra u_{\tau}')

subplot(2,2,1)
plot(k(1:round(nx/2)),power_u_phi_kx(1:round(nx/2))/max(power_u_phi_kx(1:round(nx/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('kx')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra u_{\phi}')

subplot(2,2,3)
plot(ky(1:round(ny/2)),power_u_phi_ky(1:round(ny/2))/max(power_u_phi_ky(1:round(ny/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('ky')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra u_{\phi}')

figure(90)
plot(k(1:round(nx/2)),power_phi(1:round(nx/2))/max(power_phi(1:round(nx/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('kx')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra En_{\phi}')

figure(100)
plot(ky(1:round(ny/2)),power_phi_ky(1:round(ny/2))/max(power_phi_ky(1:round(ny/2)))); hold on;
%loglog(k(1:nx/2),1e-32*k(1:nx/2).^(-3)); 
ylabel('abs(fft(rhs))^2')
xlabel('ky')
xlim([ky(1) ky(round(ny/2))])
title('1D Spectra En_{\phi}')





%figure(3)
%contourf(Lat,Lon,rhs(:,:,end)); colorbar
%xlabel('lon'); ylabel('lat')
%title('RHS')
%pbaspect([3 1 1])

figure(4)
loglog(k(1:nx/2),power_uv(1:nx/2)); 
ylabel('abs(fft(En_{uv}))^2')
xlabel('k')
xlim([k(1) k(nx/2)])



figure(5)
loglog(k(1:nx/2),power_uv(1:nx/2)); hold on;
loglog(k(1:nx/2),power_w(1:nx/2)); hold on
loglog(k(1:nx/2),power_phi(1:nx/2)); hold on
loglog(k(1:nx/2),power_tau(1:nx/2)); hold on
loglog(k(1:nx/2),1e10*k(1:nx/2).^(-3)); hold on;
%loglog(k(1:nx/2),1e35*power_rhs(1:nx/2));
legend('En_{uv500}','En_{w500}','En_{\phi}','En_{\tau}','k^{-3}','Location','NorthEast')
ylabel('Amp')
legend boxoff
xlabel('k')
xlim([k(1) k(nx/2)])
title('1-D Spectra GCM')

