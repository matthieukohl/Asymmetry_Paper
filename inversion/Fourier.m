% Calculates the Power-Spectrum of RHS and Omega_QG

clear;
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;

r=1;

for i=3:16

path_temp = strcat('/disk7/mkohl/interpolation/output1/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/u_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/v_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');


lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
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
event_lonspan = [1:256]; % 1:256
event_timespan = [1:100];
level_indices = [5:30]; % 5-30;
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
u_in= ncread(path_u,'temp_plevel', start, count);
v_in= ncread(path_v,'temp_plevel', start, count);
omega_in= ncread(path_omega,'temp_plevel', start, count);

[T, u, v, omega] = ...
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
            
[A, B, C, J, J1, J2, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

              
% spatially varying static stability: sigma_accu
% this quantity is mainly used to calculate J term
Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;


% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);


% term A (Q vector)

for k = 1 : length(level)

    % probably the ug and vg field need smoothing as well...
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

for k = 1:length(level)
    
    RHS(:,k,t) = fftshift(fft(rhs(15,:,k,t)));
    RHS_power_place(:,k,t) = abs(RHS(:,k,t)).^2;
    V(:,k,t) = fftshift(fft(v(15,:,k,t)));
    V_power_place(:,k,t) = abs(V(:,k,t)).^2;
    U(:,k,t) = fftshift(fft(u(15,:,k,t)));
    U_power_place(:,k,t) = abs(U(:,k,t)).^2;
    W(:,k,t) = fftshift(fft(omega(15,:,k,t)));
    W_power_place(:,k,t) = abs(W(:,k,t)).^2;
end

RHS_power_place = squeeze(mean(RHS_power_place(:,:,:),2));
V_power_place = squeeze(mean(V_power_place(:,:,:),2));
U_power_place = squeeze(mean(U_power_place(:,:,:),2));
W_power_place = squeeze(mean(W_power_place(:,:,:),2));

end

RHS_power(:,:,i-2) = RHS_power_place;
V_power(:,:,i-2) = V_power_place;
U_power(:,:,i-2) = U_power_place;
W_power(:,:,i-2) = W_power_place;

end

%% Plot Fourier Spectrum 

RHS_power = reshape(RHS_power,[length(lon),length(event_timespan)*14]);
V_power = reshape(V_power,[length(lon),length(event_timespan)*14]);
U_power = reshape(U_power,[length(lon),length(event_timespan)*14]);
W_power = reshape(W_power,[length(lon),length(event_timespan)*14]);

RHS_power = mean(RHS_power(:,:),2);
V_power = mean(V_power(:,:),2);
U_power = mean(U_power(:,:),2);
W_power = mean(W_power(:,:),2);

R_eff = cosd(lat(15))*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);
Wavenumber = R_eff*k;

Theoretical = 5*10^(6)*Wavenumber.^(-3);

k_centroid_V = sum(k(129:214).*V_power(129:214)')/sum(V_power(129:214));
k_centroid_RHS = sum(k(129:214).*RHS_power(129:214)')/sum(RHS_power(129:214));
k_centroid_W = sum(k(129:214).*W_power(129:214)')/sum(W_power(129:214));

n_centroid_V = k_centroid_V*R_eff;
n_centroid_RHS = k_centroid_RHS*R_eff;
n_centroid_W = k_centroid_W*R_eff;


figure(1)

plot(Wavenumber(129:214),RHS_power(129:214));
xlabel('Wavenumber (n)')
ylabel('|RHS|^2')
title('Power Spectrum RHS at lat = 45\circ r=0.1')

saveas(figure(5),'PowerW.png');

figure(2)

loglog(Wavenumber(129:214),V_power(129:214)); hold on;
loglog(Wavenumber(136:163), Theoretical(136:163));
xlabel('Wavenumber (n)')
ylabel('|V|^2')
title('Power Spectrum V at lat = 45\circ r=0.1')

% figure(3)
% 
% loglog(Wavenumber(129:218),U_power(129:218)); hold on;
% loglog(Wavenumber(136:163), Theoretical(136:163));
% xlabel('Wavenumber (n)')
% ylabel('|U|^2')
% title('Power Spectrum U at lat = 45\circ r=0.1')

% figure(4)
% 
% loglog(Wavenumber(129:218),U_power(129:218)+V_power(129:218)); hold on;
% loglog(Wavenumber(136:163), Theoretical(136:163));
% xlabel('Wavenumber (n)')
% ylabel('|U|^2+|V|^2')
% title('Power Spectrum U^2+V^2 at lat = 45\circ r=0.1')

figure(5)

plot(Wavenumber(129:214),W_power(129:214)); hold on;
%loglog(Wavenumber(136:163), Theoretical(136:163));
xlabel('Wavenumber (n)')
ylabel('|W|^2')
title('Power Spectrum W^2 at lat = 45\circ r=0.1')

