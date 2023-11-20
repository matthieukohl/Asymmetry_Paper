clear;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/temp_day0075h00.nc';
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;

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

omega_QG_RHS = [];
omega_QG_LHS = [];
omega = [];

    


for i=1:16
    
    s = load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    %t = load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/LHS_rhs_smooth0_iter10_day',num2str(25 * i),'.mat')); 
    omega_QG_RHS=cat(4,omega_QG_RHS,s.omega_QG);
    %omega_QG_LHS = cat(4,omega_QG_LHS,t.omega_QG);
end
 
 
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);
skewness_omegaQG_RHS = Skewness(omega_QG_RHS);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS);



%% Total Contribution to Skewness

omega_QG_total= [];
sigma_new= [];
sigma_new_smoothed= [];

for i=1:16
    
    t= load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/random/varyingf/complete/Total_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    omega_QG_total=cat(4,omega_QG_total,t.omega_QG);
    omega = cat(4,omega,t.omega);     
end


omega_QG_total = omega_QG_total-mean(omega_QG_total(:,:,:,:),2);
omega = omega-mean(omega(:,:,:,:),2);


skewness_omega = Skewness(omega);
lambda_omega = Lambda(omega);

skewness_omegaQG_total = Skewness(omega_QG_total);
lambda_omegaQG_total = Lambda(omega_QG_total);




%% Plotting


% figure(1)
% 
% plot(skewness_omega); hold on
% %plot(skewness_omegaQG_RHS); hold on
% %plot(skewness_omegaQG_total);
% title('Skewness Vertical Velocity (r=0.1)')
% legend('omega','omegaQG_{total}');
% xlabel('time(quarter days)');
% ylabel('skewness')



figure(2)

t=10:1200;

plot(t,lambda_omega(10:1200)); hold on 
plot(t,lambda_omegaQG_total(10:1200)); hold on;
plot(t,lambda_omegaQG_RHS(10:1200)); hold on

title('Lambda Vertical Velocity (r=0.1)')
legend('omega','omegaQG_{total}','omegeaQG_{RHS}');
xlabel('time(quarter days)');
ylabel('lambda')


saveas(gcf,'Lambda01.png')




