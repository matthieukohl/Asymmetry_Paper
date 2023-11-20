% Calculates the Power-Spectrum of RHS and Omega_QG

clear;

path_temp = '/disk7/mkohl/interpolation/output_abs/abs1.0_T85/temp_day1450h00.nc';
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;
R = 6371000.0;

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);

lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;
plevels = double(ncread(path_temp, 'level'));
level_indices = [5:30];
level = plevels(level_indices);


%event_latspan = [83:111]; 
event_latspan = [93:108]; 
event_lonspan = [1:256];
event_timespan = [100:200];

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

abs_T85 = [0.4;0.6;0.7;0.8;0.9;1.0;1.2;1.4;1.6;1.8;2.0;2.5;3.0;4.0;6.0];

strat_level_turbulent = [8,9,10,10,11,11,12, 12, 13, 13, 14, 14, 15,16,18];

n_RHS = zeros(15,1);
n_W = zeros(15,1);

for jj = 1:15

jj

omega_QG = [];
rhs = [];
omega = [];
r = [];

for i= 5:6
s = load(strcat('output/abs_simulation/',list{jj},'/final/Total_rhs_smooth0_rms4_geostrophic_w_trial_day',num2str(5 *(28+i)),'.mat'));
omega= cat(4,omega,s.omega);
rhs = cat(4,rhs,s.rhs);
end

omega = omega(11:26,:,:,:);
rhs = rhs(11:26,:,:,:);

RHS = zeros(length(lon),length(level),1200);
RHS_power_place = zeros(length(lon),length(level),1200);
W = zeros(length(lon),length(level),1200);
W_power_place = zeros(length(lon),length(level),1200);

RHS = zeros(length(lat),length(lon),length(level),size(omega,4));
RHS_power_place = zeros(length(lat),length(lon),length(level),size(omega,4));
W = zeros(length(lat),length(lon),length(level),size(omega,4));
W_power_place = zeros(length(lat),length(lon),length(level),size(omega,4));

k_centroid_RHS = zeros(length(lat),1);
k_centroid_W = zeros(length(lat),1);

n_centroid_RHS = zeros(length(lat),1);
n_centroid_W = zeros(length(lat),1);

for j = 1:length(lat)
    
  for t = 1 : size(omega,4) 
    for k = 1:strat_level_turbulent(jj)
    RHS(j,:,k,t) = fftshift(fft(rhs(j,:,k,t)));
    RHS_power_place(j,:,k,t) = abs(RHS(j,:,k,t)).^2;
    W(j,:,k,t) = fftshift(fft(omega(j,:,k,t)));
    W_power_place(j,:,k,t) = abs(W(j,:,k,t)).^2;
    end
  end
    
RHS_power = squeeze(mean(mean(RHS_power_place(j,:,:,:),3),4))';
W_power = squeeze(mean(mean(W_power_place(j,:,:,:),3),4))';


%% Plot Fourier Spectrum 

R_eff = cosd(lat(j))*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);
Wavenumber = R_eff*k;

k_centroid_RHS(j) = squeeze(sum(k(129:214).*RHS_power(129:214)')/sum(RHS_power(129:214)));
k_centroid_W(j) = squeeze(sum(k(129:214).*W_power(129:214)')/sum(W_power(129:214)));

n_centroid_RHS(j) = k_centroid_RHS(j)*R_eff;
n_centroid_W(j) = k_centroid_W(j)*R_eff;

end

n_RHS(jj) = nanmean(n_centroid_RHS);
n_W(jj) = nanmean(n_centroid_W);


% h1 = figure(1);
% 
% plot(Wavenumber(129:214),RHS_power(129:214));  hold on;
% xlabel('Wavenumber (n)')
% ylabel('|RHS|^2')
% title('Power Spectrum RHS at lat = 45\circ')
% 
% %  
%h2 = figure(2);

% plot(Wavenumber(129:214),W_power(129:214)); hold on;
% xlabel('Wavenumber (n)')
% ylabel('|W|^2')
% title('Power Spectrum W^2 at lat = 45\circ abs4.0')
% savefig('Trial.fig')

end
% saveas(h1,'Spec_RHS.m');
% saveas(h2,'Spec_W.m')
% 
% figure(3)
% plot(abs_T85,n_RHS);
% xlabel('abs')
% ylabel('Wavenumber')
% title('Centroid Wavenumber RHS-Specrum')
% saveas(gcf,'RHS.m')
% 
% 
% figure(4)
% plot(abs_T85,n_W)
% xlabel('abs')
% ylabel('Wavenumber')
% title('Centroid Wavenumber W-Specrum')
% saveas(gcf,'W.m')
% 







