% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS from reanalysis

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

tic

% Define Path 

% high resolution

path = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_atlantic_patch_inversion_with_uv.nc';

% coarse resolution: this uses the ERA5 2.5 fields that are not truly
% coarse grained

%path = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_atlantic_patch_inversion_with_uv_coarse_resolution_1p4.nc';

% choose if coarse grain or not

resolution = 'high';
tf_resolution = strcmp(resolution,'coarse');

% Read in Grid-Variables Lat/Lon 

lat = ncread(path,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(path,'longitude');lon = double(lon);

level = ncread(path,'level'); level = double(level)*100; level = flip(level);
%level = level(1:30);
time = ncread(path,'time');

% Read in Variables (Lat/Lon need to be flipped)

start = [1, 1, 1, 1];
count = [length(lon), length(lat), length(level), length(time)];

T_in= ncread(path,'t', start, count);
u_in= ncread(path,'u', start, count);
v_in= ncread(path,'v', start, count);
omega_in= ncread(path,'w', start, count);

[u,v,T,omega] = deal(zeros(length(lat), length(lon), length(level), length(time)));

for t = 1:length(time)
u(:,:,:,t) = rearrange_era(u_in(:,:,:,t),lat,lon,level);
v(:,:,:,t) = rearrange_era(v_in(:,:,:,t),lat,lon,level);
T(:,:,:,t) = rearrange_era(T_in(:,:,:,t),lat,lon,level);
omega(:,:,:,t) = rearrange_era(omega_in(:,:,:,t),lat,lon,level);
end

clear('u_in','v_in','omega_in','T_in');

% option to coarse grain

res = 10;

if tf_resolution==1
% average
u = movmean(movmean(u,res+1,1),res+1,2);
v = movmean(movmean(v,res+1,1),res+1,2);
T = movmean(movmean(T,res+1,1),res+1,2);
omega = movmean(movmean(omega,res+1,1),res+1,2);
% downsample
u = u(1:res:end,1:res:end,:,:);
v = v(1:res:end,1:res:end,:,:);
T = T(1:res:end,1:res:end,:,:);
omega = omega(1:res:end,1:res:end,:,:);
end

if tf_resolution==1
    lat = lat(1:res:end);
    lon = lon(1:res:end);
end

dphi = abs(lat(2) - lat(1)) / 180.0 * 3.1415926;
dlambda = abs(lon(2) - lon(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

event_latspan = [1:length(lat)]; %83 -111
event_lonspan = [1:length(lon)]; % 1:256
event_timespan = [1:length(time)];
level_indices = [1:length(level)];


phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% pre-allocation

[A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

err_rms = zeros(length(event_timespan),1);


% varying f
       
% f = 2 * Omega * sind(lat);
% beta = 1 / R * d_dphi(f, dphi);

% constant f0
lat0 = lat(1)+(lat(end)-lat(1))/2;
f = 2 * Omega * sind(lat0*ones(size(lat)));
beta = 2 * Omega * cosd(lat0*ones(size(lat)))/R;
    
for t = 1 : length(event_timespan)

t    
    
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r_pholder = r_parameter1(T(:,:,:,t),Level,level);

% rescale r by 1/10 to test how high lambda gets for reanalysis

%r_pholder(r_pholder<1) = 0.1*r_pholder(r_pholder<1);

r(:,:,:,t) = r_pholder;

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
    
while rms(diff(:))>=1*10^(-4) && count<400
    
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
%error(count) = rms(diff(:));
% asym(count) = Lambda_weighted(w_new(11:26,:,1:strat_level_turbulent(tt),1),lat(11:26));
% 
% plot(squeeze(w_new(15,:,10))); hold on;
                                                                                                                  
end
count

err_rms(t) = rms(diff(:));
omega_QG(:,:,bound:end,t) = w_new;
end

toc

% figure(6)
% plot(lon,omega(11,:,7+bound,1)); hold on;
% plot(lon,w_new(11,:,7,1)); hold on;
% plot(lon,test(11,:,7,1)); hold on;
% legend('omega','QG','dry')
% xlabel('lon')
% ylabel('w(Pa/s)')
% title('w at 500hPA and 40 lat, r=0.4')

if tf_resolution==1
    
save(strcat(['/net/graupel/archive1/users/mkohl/inversion/output/reanalysis/rms1em4_geostrophic_boundary_w_without_rescale_coarse_resolution_res',num2str(res),'.mat']),'omega','omega_QG','r','rhs','lat','err_rms');

else
    
save(strcat('/net/graupel/archive1/users/mkohl/inversion/output/reanalysis/rms1em4_geostrophic_boundary_w_without_rescale.mat'),'omega','omega_QG','r','rhs','lat','err_rms');

end


[a,p500] = min(abs(level-500*1e2));
[k_500_new,spec_500_true] = Spectrum2(omega(:,:,p500,t),lat,lon);
[k_500_new,spec_500_inverted] = Spectrum2(omega_QG(:,:,p500,t),lat,lon);

R_eff = R*cosd(45);
N = length(lon);
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
Ld = 1e6/(2*sqrt(2));
k = k*Ld;

figure
loglog(k(1:N/2),spec_500_true(1:N/2)/max(spec_500_true(1:N/2))); hold on;
loglog(k(1:N/2),spec_500_inverted(1:N/2)/max(spec_500_inverted(1:N/2))); hold on;
loglog(k(1:N/2),1e1*k(1:N/2).^(-3),'k--'); hold on;
xlabel('Wavenumber k');
ylabel('Spec')
legend('true','inverted','k^{-3}','Location','SouthWest'); legend boxoff


% plot the fields

cmax = 1;
nint = 20;
cint = cmax/nint;

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

subplot(2,1,1)

contourf(-omega(:,:,p500,t),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega')

subplot(2,1,2)

contourf(-omega_QG(:,:,p500,t),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega Inverted')

% calculate lambda

% lambda_true = Lambda_weighted(omega(:,:,p500,t)-mean(mean(omega(:,:,p500,t))),lat)
% 
% lambda_inverted = Lambda_weighted(omega_QG(:,:,p500,t)-mean(mean(omega_QG(:,:,p500,t))),lat)

lambda_true = Lambda_weighted(omega(:,:,p500,t),lat)

lambda_inverted = Lambda_weighted(omega_QG(:,:,p500,t),lat)


lambda_cal = [lambda_true,lambda_inverted];

bar(lambda_cal)
labels={'true'; 'inverted'};
set(gca,'xticklabel',labels)
ylabel('$\lambda$')
set(gca,'FontSize',12)