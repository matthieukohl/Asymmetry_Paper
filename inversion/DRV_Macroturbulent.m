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

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

%list = {'abs3.0_T85','abs4.0_T85','abs6.0_T85'};

%list = {'abs1.8_T85','abs2.0_T85','abs2.5_T85'};

%list = {'abs1.2_T85','abs1.4_T85','abs1.6_T85'};

%list = {'abs0.8_T85','abs0.9_T85','abs1.0_T85'};

%list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85'};

strat_level_turbulent = [8,9,10,10,11,11,12, 12, 13, 13, 14, 14, 15,16,18];

for tt = 6:6
tt
for i= 4:4
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/temp_day',num2str(5 *(28+i)),'0h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/ug_day',num2str(5 *(28+i)),'0h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/vg_day',num2str(5 *(28+i)),'0h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output_abs/',list{tt},'/omega_day',num2str(5 *(28+i)),'0h00.nc');

%path_temp = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/temp_day0560h00.nc');
%path_u = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/ug_day0560h00.nc');
%path_v = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/vg_day0560h00.nc');
%path_omega = strcat('/disk7/mkohl/interpolation/output_abs_mode/',list{tt},'/omega_day0560h00.nc');


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

event_latspan = [83:113]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [1:200];
%event_timespan = [6:9];
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
       
f = 2 * Omega * sind(lat);
beta = 1 / R * d_dphi(f, dphi);

% f0 = 2 * Omega * sin(lat_series_original(100) / 180 * 3.1415926);
% beta = 2 * Omega / R * cos(lat_series_original(100) / 180 * 3.1415926);
            
[A, B, C, zeta, J, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[r,density,lapse] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
[omega_QG,omega_QG_RHS] = deal(nan(length(event_latspan),length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

% static stability: sigma
% the following formula is more accurate since it takes into account the 
% moist part in the exponent when calculating the potential temperature

Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level).^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

% T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));
% theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
% sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);

r(:,:,:,t) = r_parameter1(T(:,:,:,t),Level,level);

density(:,:,:,t) = rho(T(:,:,:,t),Level);

lapse(:,:,:,t) = -9.81*density(:,:,:,t).*d_dp(T(:,:,:,t),level)*10^3; % K/km

for k = 1:length(level)
    vor(:,:,k,t) = curl(phi, u(:,:,k,t), v(:,:,k,t), dphi, dlambda);
end

end


end

end

% dp = level(6)-level(16);
% p = level(6:16);
% N2 = squeeze(sigma_accu(27,77,6:16,4));
% rz = squeeze(r(27,77,6:16,4));
% S0_eff = trapz(p,2/(pi*dp)*sin(pi*(level(6)-p)/dp).*N2.*rz);
% S1_eff = trapz(p,1/(2*dp)*(1-4/pi*sin(pi*(level(6)-p)/dp)).*N2.*rz);
% S2_eff = trapz(p,1/(2*dp)*cos(pi*(level(6)-p)/dp).*N2.*rz);
% S0 =  trapz(p,2/(pi*dp)*sin(pi*(level(6)-p)/dp).*N2);
% r_eff = S0_eff/S0;
% 
% figure(1); plot(squeeze(omega(27,77,:,4)),level); set(gca,'Ydir','reverse')
% 
% figure(2); plot(squeeze(r(27,77,:,4)),level); set(gca,'Ydir','reverse')


% % Find Paul's DRV: lat = 61, lon = 158, i = 4; t = 139
% [Lat,Lon] = meshgrid(lat,lon);
% Lat = Lat'; Lon = Lon';
% contourf(Lon,Lat,omega(:,:,12,139));colorbar;

% find maximum -omega

time = 139;

[latm,lonm] = find(omega(:,:,12,time)==min(min(omega(:,:,12,time))));

% find troposphere

lapse_mean = mean(mean(mean(lapse,1),2),4);

trop = find(lapse_mean>-2,1);

for jj = 1:length(lon)


% find lowest level where r>0

surf = find(r(latm,jj,:,time)>0,1);

% Test Theories

dp = level(surf)-level(trop);
p = level(surf:trop);
N2 = squeeze(sigma_accu(latm,jj,surf:trop,time));
rz = squeeze(r(latm,jj,surf:trop,time));
S0_eff(jj) = trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2.*rz);
S0(jj) =  trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2);
r_eff(jj) = S0_eff(jj)/S0(jj);

r_mean(jj) = mean(rz);

end

%%
figure(10);
plot(lon,r_eff); hold on;
plot(lon,-omega(latm,:,12,time)); 
xlim([lon(1) lon(end)])
xlabel('lon'); legend('r_{eff}','-\omega_{430}')
legend boxoff
title('Macroturbulent abs4.0')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'fontsize', 14);

% find lowest level where r>0

surf = find(r(latm,lonm,:,1)>0,1);

% Test Theories

dp = level(surf)-level(trop);
p = level(surf:trop);
N2 = squeeze(sigma_accu(latm,lonm,surf:trop,time));
rz = squeeze(r(latm,lonm,surf:trop,time));
S0_eff = trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2.*rz);
S1_eff = trapz(p,1/(2*dp)*(1-4/pi*sin(pi*(level(surf)-p)/dp)).*N2.*rz);
S2_eff = trapz(p,1/(2*dp)*cos(pi*(level(surf)-p)/dp).*N2.*rz);
S0 =  trapz(p,2/(pi*dp)*sin(pi*(level(surf)-p)/dp).*N2);
r_eff = S0_eff/S0;

figure(1); plot(squeeze(omega(latm,lonm,:,time)),level/10^2,'linewidth',1.5); 
set(gca,'Ydir','reverse'); xlabel('omega'),ylabel('p (hPA)');
ylim([level(trop)/10^2 level(surf)/10^2])
title('Macroturbulent')
set(gca,'fontsize', 14); 


figure(2); plot(squeeze(r(latm,lonm,:,time)),level/10^2,'linewidth',1.5); 
set(gca,'Ydir','reverse')
xlabel('r'),ylabel('p (hPA)'); ylim([level(trop)/10^2 level(surf)/10^2])
title('Macroturbulent')
set(gca,'fontsize', 14); 







% 
% 
% plot(squeeze(r(26,158,:,139)),level/10^2); set (gca,'Ydir','reverse'); 
% xlabel('Reduction factor r'); ylabel('p(hPA)')

%r1 = min(r,0);
%r1(r1==0) = NaN;



%levs = linspace(min(min(r(:,:,16,130))),-min(min(r(:,:,16,139))),20);

% %Mask All r values where omega>0
% 
% for ii = 1:length(lat)
%     for jj = 1:length(lon)
%         for kk = 1:length(level)
%             for tt = 1:length(event_timespan)
%         
%         if omega(ii,jj,kk,tt)>0
%             r(ii,jj,kk,tt) = NaN;
%         else
%         end
%             end
%         end
%     end
% end
% 
% 
% 
% figure('pos', [10, 10, 1000, 300]) 
% contourf(Lon,Lat,r(:,:,16,139),'LineColor', 'none'); colorbar
% caxis([min(min(r(:,:,16,139))) -min(min(r(:,:,16,139)))])
% colormap(redblue)
% set(gca,'Color',[0.5 0.5 0.5])
% xlabel('lon'); ylabel('lat'); title(['Reduction factor r at  ', num2str(round(level(16)/10^2)),'hPA'])
% 
% 
% figure('pos', [10, 10, 1000, 300]) 
% contourf(Lon,Lat,r(:,:,11,139),'LineColor', 'none'); colorbar
% caxis([min(min(r(:,:,11,139))) -min(min(r(:,:,11,139)))])
% colormap(redblue)
% set(gca,'Color',[0.5 0.5 0.5])
% xlabel('lon'); ylabel('lat'); title(['Reduction factor r at  ', num2str(round(level(11)/10^2)),'hPA'])
% % 
% figure('pos', [10, 10, 1000, 300]) 
% contourf(Lon,Lat,r(:,:,19,139),'LineColor', 'none'); colorbar
% caxis([min(min(r(:,:,19,139))) -min(min(r(:,:,19,139)))])
% colormap(redblue)
% set(gca,'Color',[0.5 0.5 0.5])
% xlabel('lon'); ylabel('lat'); title(['Reduction factor r at  ', num2str(round(level(19)/10^2)),'hPA'])



% Mode

% figure('pos', [10, 10, 1000, 300]) 
% contourf(Lon,Lat,r(:,:,16,end),'LineColor', 'none'); colorbar
% caxis([min(min(r(:,:,16,end))) -min(min(r(:,:,16,end)))])
% colormap(redblue)
% xlabel('lon'); ylabel('lat'); title(['Reduction factor r at  ', num2str(round(level(16)/10^2)),'hPA'])




%figure(2); plot(squeeze(r(12,142,:,4)),level/10^2,'LineWidth',1.6); 
%xlabel('Reduction factor r'); ylabel('p(hPA)'); set (gca,'Ydir','reverse')
%saveas(gcf,'/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/r_DRV_mode','epsc');


% v = VideoWriter('/disk7/mkohl/Mode_Zurita_Gotor/Plots_Proposal_Main/Large_Domain/DRV_Macroturbulence.avi');
% v.Quality = 100;
% v.FrameRate = 8;
% open(v)
% figure('pos', [10, 10, 1000, 300])
% for tt = 127:140
%     contourf(Lon(:,90:130),Lat(:,90:130),omega(:,90:130,10,tt)); colorbar; % caxis([-1.6 0.2])
%     xlabel('lon'); ylabel('lat')
%     frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% close(v)

