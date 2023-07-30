% Calculates the Power-Spectrum of RHS and Omega_QG

clear;
close all;

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale1.0/temp_day0125h00.nc';
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


event_latspan = [83:111]; 
event_lonspan = [1:256];
event_timespan = [1:100];

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};

r = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

%strat_level = [19;14;13;14;14;13;12;18;12;18];

strat_level = [18,18,18,18,17,17,16,16,16,16]-4;

n_RHS = zeros(10,1);
n_W = zeros(10,1);

for jj = 1:10
    
jj
    
omega_QG = [];
rhs = [];
omega = [];

for ii= 1:4
s = load(strcat('output/r_simulation/rescale',list{jj},'/final/Total_rhs_smooth0_rms4_geostrophic_day',num2str(25 *(8+ii)),'.mat'));
omega= cat(4,omega,s.omega);
rhs = cat(4,rhs,s.rhs);
end

omega = omega(:,:,5:30,:);
rhs = rhs(:,:,5:30,:);

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
    for k = 1:strat_level(jj)
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

end
% saveas(h1,'Spec_RHS_r.m');
% saveas(h2,'Spec_W_r.m')
% 
% % 

% figure(3)
% plot(r,n_RHS);
% xlabel('abs')
% ylabel('Wavenumber')
% title('Centroid Wavenumber RHS-Specrum all lat')
% savefig('RHS_r_full.fig')
% 
% 
% figure(4)
% plot(r,n_W)
% xlabel('abs')
% ylabel('Wavenumber')
% title('Centroid Wavenumber W-Specrum all lat')
% savefig('W_r_full.fig')

%Wavenumber = [18.5712, 13.7374, 10.6823, 14.2371, 9.8514, 6.1367, 8.9836,....
%8.0164, 6.7059, 7.9694];