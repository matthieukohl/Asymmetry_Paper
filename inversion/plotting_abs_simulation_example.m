% plotting abs simulation example

% modal regime

clear;
close all;

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

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};
tic
for tt = 6:6
    
tt    

% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/temp_day0560h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/ug_day0560h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/vg_day0560h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/omega_day0560h00.nc');

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
event_timespan = [9:9]; % 6-9
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
[r,rhs, omega_b] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)
    
% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
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

w_old = randn(length(event_latspan), length(event_lonspan), length(level_place)); 

% w_new = SIP_diagnostic1(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
% sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
sigma_place(:), dphi, dlambda, rhs_place(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));

diff = w_new - w_old;

% Pre-Allocation

sigma_r = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));
del2_J = deal(zeros(length(event_latspan),length(event_lonspan),length(level_place)));

%while rms(diff(:))>=10^(-15) && count<350
while rms(diff(:))/rms(w_new(:))>=10^(-3) && count<350
rms(diff(:));
count = count + 1;
       
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
%                    sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));

w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));                            
                             
                             
diff = w_new - w_old; 
error(count) = rms(diff(:));
asym(count) = Lambda(w_new(:,:,1:8));
end

omega_QG(:,:,bound:end,t) = w_new;

end

end

omega_mode = omega(:,:,:,t);
omega_qg_mode = omega_QG(:,:,:,t);

% macroturbulent runs

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};


for tt = 6:6
tt
for i= 6:6 %1:1 %2:2 20/20
    i
err_rms = zeros(200,1);
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/ug_day',num2str(5 *(28+i)),'0h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/vg_day',num2str(5 *(28+i)),'0h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');

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
event_timespan = [150:150];
%event_timespan = [10:10];
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
%u_in= ncread(path_u,'u_plevel', start, count);
%v_in= ncread(path_v,'v_plevel', start, count);
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

%w_new = randn(size(w_old));

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
    
%while rms(diff(:))>=1*10^(-4) && count<50
while rms(diff(:))/rms(w_new(:))>=10^(-3) && count<350    

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

end

end

omega_turb = omega(:,:,:,t);
omega_qg_turb = omega_QG(:,:,:,t);

[Lon,Lat] = meshgrid(lon,lat);
[a,ind500] = min(abs(level-50000));

fig = figure('pos', [10, 10, 1200, 500]);

set(groot,'DefaultTextInterpreter','default');
set(groot,'DefaultLegendInterpreter','default');

nlevels = 11;
cint = max(max(abs(omega_mode(:,:,ind500))))/nlevels;

s1 = subplot(2,2,1);
contourf(lon,lat,-omega_mode(:,:,ind500),[-cint*nlevels:cint:cint*nlevels],'EdgeColor','none'); 
%colormap(redbluecmap(nlevels));
%hc=colorbar;
%title(hc,'Pa/s');
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel(sprintf('Latitude (%c)', char(176)))
%t1 = title('\rm GCM Vertical Velocity');
%t1pos = get(t1,'position');
%set(t1,'position',t1pos+[0 1 0]);
set(gca, 'TickDir', 'out')
set(gca,'linewidth',1.5)
caxis([-cint*nlevels cint*nlevels])

s2 = subplot(2,2,2);
contourf(lon,lat,-omega_qg_mode(:,:,ind500),[-cint*nlevels:cint:cint*nlevels],'EdgeColor','none'); 
%colormap(redbluecmap(nlevels));
hc=colorbar;
hcpos = get(hc,'Position');
set(hc,'Units','normalized', 'position', hcpos + [0.01 0 0 0 ]);
%cbpos = get(hc,'Position');
%cbpos(3) = 0.4*cbpos(3); % narrower
%cbpos(1) = cbpos(1)+0.09; 
%a = hc.Position;
%set(hc,'Position',[a(1)+0.1 a(2) a(3) a(4)])
title(hc,'Pa s^{-1}');
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel(sprintf('Latitude (%c)', char(176)))
%t1 = title('\rm QG Vertical Velocity');
%t1pos = get(t1,'position');
%set(t1,'position',t1pos+[0 0.04 0]);
set(gca, 'TickDir', 'out')
set(gca,'linewidth',1.5)
caxis([-cint*nlevels cint*nlevels])

s1Pos = get(s1,'position');
s2Pos = get(s2,'position');
s2Pos(3:4) = [s1Pos(3:4)];
set(s2,'position',s2Pos);

set(s1,'Units','normalized', 'position', s1Pos - [0.05 0 0 0]);
set(s2,'Units','normalized', 'position', s2Pos - [0.05 0 0 0]);

% cm = redbluecmap; %cm(6,:) = [1 1 1]; colormap(cm);
% cm = imresize(cm, [nlevels, 3]); cm = min(max(cm, 0), 1);
% cm(0.5*(nlevels+1),:) = [1 1 1]; 
% colormap(cm)
% set(fig,'renderer','Painters')

nlevels = 11;
cint = max(max(abs(omega_turb(:,:,ind500))))/nlevels;


s3 = subplot(2,2,3);
contourf(lon,lat,-omega_turb(:,:,ind500),[-cint*nlevels:cint:cint*nlevels],'EdgeColor','none'); 
%contourf(lon,lat,-omega_turb(:,:,ind500),[-1:0.05:1],'EdgeColor','none'); 


%colormap(redblue(24));
%hc=colorbar;
%title(hc,'Pa/s');
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel(sprintf('Latitude (%c)', char(176)))
%title('\rm QG')
set(gca, 'TickDir', 'out')
set(gca,'linewidth',1.5)
caxis([-cint*nlevels cint*nlevels])
%caxis([-1 1])

s4 = subplot(2,2,4);
%cmlev = 15;
contourf(lon,lat,-omega_qg_turb(:,:,ind500),[-cint*nlevels:cint:cint*nlevels],'EdgeColor','none'); 
%contourf(lon,lat,-omega_qg_turb(:,:,ind500),[-1:0.05:1],'EdgeColor','none'); 
%redbluecmap(nlevels)
%cm = redbluecmap(nlevels);
%cm = imresize(cm, [cmlev, 3]); cm = min(max(cm, 0), 1);
%cm(0.5*(cmlev+1),:) = [1 1 1]; 
%colormap(cm)
%cm(6,:) = [1 1 1];
%colormap(redblue(24))
hc=colorbar;
title(hc,'Pa s^{-1}');
hcpos = get(hc,'Position');
set(hc,'Units','normalized', 'position', hcpos + [0.01 0 0 0 ]);
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel(sprintf('Latitude (%c)', char(176)))
%title('\rm QG')
set(gca, 'TickDir', 'out')
set(gca,'linewidth',1.5)
caxis([-cint*nlevels cint*nlevels])
%caxis([-1 1])

s3Pos = get(s3,'position');
s4Pos = get(s4,'position');
s4Pos(3:4) = [s3Pos(3:4)];
set(s4,'position',s4Pos);

set(s3,'Units','normalized', 'position', s3Pos - [0.05 0 0 0]);
set(s4,'Units','normalized', 'position', s4Pos - [0.05 0 0 0]);

% a,b,c,d
annotation('textbox', [0.04, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none')
annotation('textbox', [0.48, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none')
annotation('textbox', [0.04, 0.47, 0.01, 0.03], 'String', "(c)",'FontSize',12,'EdgeColor','none')
annotation('textbox', [0.48, 0.47, 0.01, 0.03], 'String', "(d)",'FontSize',12,'EdgeColor','none')

% text boxes

annotation('textbox', [0.17, 0.91, 0.3, 0.1], 'String', "GCM Vertical Velocity",'EdgeColor','none','FontSize',12)
annotation('textbox', [0.61, 0.91, 0.3, 0.1], 'String', "Inverted Vertical Velocity",'EdgeColor','none','FontSize',12)

annotation('textbox', [0.91, 0.8, 0.01, 0.01], 'String', "Unstable Mode",'EdgeColor','none','FontSize',12)
annotation('textbox', [0.91, 0.3, 0.01, 0.01], 'String', "Macro- Turbulence",'EdgeColor','none','FontSize',12)


% with two colorbars

% % a,b,c,d
% annotation('textbox', [0.07, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.52, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.07, 0.47, 0.01, 0.03], 'String', "(c)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.52, 0.47, 0.01, 0.03], 'String', "(d)",'FontSize',12,'EdgeColor','none')
% 
% % text boxes
% 
% annotation('textbox', [0.21, 0.91, 0.3, 0.1], 'String', "GCM Vertical Velocity",'EdgeColor','none','FontSize',12)
% annotation('textbox', [0.65, 0.91, 0.3, 0.1], 'String', "Inverted Vertical Velocity",'EdgeColor','none','FontSize',12)
% 
% annotation('textbox', [0.91, 0.8, 0.01, 0.01], 'String', "Unstable Modes",'EdgeColor','none','FontSize',12)
% annotation('textbox', [0.91, 0.3, 0.01, 0.01], 'String', "Macro- Turbulence",'EdgeColor','none','FontSize',12)

%colormap(redblue(nlevels))
%colormap(redbluecmap)
cm = redbluecmap; %cm(6,:) = [1 1 1]; colormap(cm);
cm = imresize(cm, [nlevels, 3]); cm = min(max(cm, 0), 1);
cm(0.5*(nlevels+1),:) = [1 1 1]; 
colormap(cm)
set(fig,'renderer','Painters')
saveas(gcf,'figures_paper/w_snapshot_2d','epsc');


% plot a slice
 
figure('pos', [10, 10, 1000, 300])
[a,lat50] = min(abs(lat-50));
 
subplot(1,2,1)
plot(lon,-omega_mode(lat50,:,ind500),'k'); hold on;
plot(lon,-omega_qg_mode(lat50,:,ind500),'b')
legend('GCM','Inversion','Location','Northwest')
legend boxoff
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel('-\omega (Pa s^{-1})')
title('\rm Modes')
xlim([lon(1) lon(end)])
set(gca,'box','off')
set(gca,'FontSize',12)


subplot(1,2,2)
plot(lon,-omega_turb(lat50,:,ind500),'k'); hold on;
plot(lon,-omega_qg_turb(lat50,:,ind500),'b')
legend('GCM','Inversion','Location','Northwest')
legend boxoff
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel('-\omega (Pa s^{-1})')
title('\rm Macroturbulence')
xlim([lon(1) lon(end)])
set(gca,'box','off')
set(gca,'FontSize',12)

annotation('textbox', [0.07, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none')
annotation('textbox', [0.51, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none')

set(gca,'FontSize',12)

subplot(1,2,2)
plot(lon,-omega_turb(lat50,:,ind500),'k'); hold on;
plot(lon,-omega_qg_turb(lat50,:,ind500),'b')
legend('GCM','Inversion','Location','Northwest')
legend boxoff
xlabel(sprintf('Longitude (%c)', char(176)))
ylabel('-\omega (Pa s^{-1})')
title('\rm Macroturbulence')
xlim([lon(1) lon(end)])
set(gca,'box','off')
set(gca,'FontSize',12)

annotation('textbox', [0.07, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none')
annotation('textbox', [0.51, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none')

saveas(gcf,'figures_paper/w_snapshot_1d','epsc')


% check altitude of max omega
% 
% omega_up = zeros(size(omega));
% 
% for ii = 1:29
%     for jj = 1:256
%         for tt = 1:150
%             if -omega(ii,jj,12,tt)>0
%                 omega_up(ii,jj,:,tt) = omega(ii,jj,:,tt);
%             else
%                 omega_up(ii,jj,:,tt) = NaN*ones(length(level),1);
%             end
%         end
%     end
% end
% 
% plot(squeeze(nanmean(nanmean(nanmean(-omega_up(:,:,:,:),1),2),4)),level/100);
% set(gca,'Ydir','reverse')
% xlabel('-\omega')
% ylabel('p(hPa)')