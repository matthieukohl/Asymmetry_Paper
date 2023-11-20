clear;

%% Concatenate Files

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0075h00.nc';
lat_series_original = ncread(path_temp, 'latitude');
event_latspan = [83:111];
lat = lat_series_original(event_latspan);

omega=[];
omega_QG=[];

for i=1:12
    
    s = load(strcat('output/rescale_1.0/sigma/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));

    omega = cat(4,omega,s.omega);
    omega_QG = cat(4,omega_QG,s.omega_QG);
    
end

omega = omega - mean(omega(:,:,:,:),2);

omega_QG = omega_QG - mean(omega_QG(:,:,:,:),2);


skewness_omega = Skewness(omega);
skewness_omegaQG = Skewness(omega_QG);

lambda_omega = Lambda(omega);
lambda_omegaQG = Lambda(omega_QG);


%% Plotting


% figure(1)
% 
% plot(skewness_omega); hold on
% plot(skewness_omegaQG); hold on
% title('Skewness Vertical Velocity (r=1)')
% legend('omega','omega_{QG}');
% xlabel('time(quarter days)');
% ylabel('skewness')
% 
% savefig('output/rescale_1.0/skewness.fig');

figure(2)

t=300:1200;

plot(t,lambda_omega(300:end)); hold on 
plot(t,lambda_omegaQG(300:end)); hold on
title('Lambda Vertical Velocity (r=1)')
legend('omega','omega_{QG}');
xlabel('time(quarter days)');
ylabel('lambda')

saveas(gcf,'Lambda10.png')

%savefig('output/rescale_1.0/lambda.fig');

%% Averages

% m = squeeze(mean(mean(omega_QG(:,:,:,:),1),2));
% n = squeeze(mean(mean(omega(:,:,:,:),1),2));
% 
% m = mean(m(:,:),1);
% n = mean(n(:,:),1);
% 
% figure(3); 
% 
% plot(m); hold on; 
% plot(n);
% title('Averages')
% legend('QG','Omega')

