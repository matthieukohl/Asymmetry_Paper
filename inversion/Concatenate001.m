clear;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0075h00.nc';
lat_series_original = ncread(path_temp, 'latitude');
event_latspan = [83:111];
lat = lat_series_original(event_latspan);

omega_QG_RHS = [];
omega = [];

for i=1:12
    
    s = load(strcat('output/rescale_0.01/sigma/equivalent_J/zero_BC/RHS_rhs_smooth0_iter10_day',num2str(25 * i),'.mat'));
    omega_QG_RHS = cat(4, omega_QG_RHS, s.omega_QG);
    omega = cat(4,omega,s.omega);   
end
    

omega = omega-mean(omega(:,:,:,:),2);
 
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);


skewness_omega = Skewness(omega);
skewness_omegaQG_RHS = Skewness(omega_QG_RHS);

lambda_omega = Lambda(omega);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS);


%% Total Contribution to Skewness

omega_QG_total=[];

for i=1:12
    
    t= load(strcat('output/rescale_0.01/sigma/equivalent_J/zero_BC/Total_rhs_smooth0_iter10_day',num2str(25 * i),'.mat'));
    
    omega_QG_total=cat(4,omega_QG_total,t.omega_QG);
    
end




omega_QG_total = omega_QG_total-mean(omega_QG_total(:,:,:,:),2);


skewness_omegaQG_total = Skewness(omega_QG_total);
lambda_omegaQG_total = Lambda(omega_QG_total);




%% Plotting


figure(1)

plot(skewness_omega); hold on
plot(skewness_omegaQG_RHS); hold on
plot(skewness_omegaQG_total);
title('Skewness Vertical Velocity (r=0.01)')
legend('omega','omegaQG_{RHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('skewness')



figure(2)

plot(lambda_omega); hold on 
plot(lambda_omegaQG_RHS); hold on
plot(lambda_omegaQG_total);
title('Lambda Vertical Velocity (r=0.01)')
legend('omega','omegaQG_{RHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('lambda')



% figure(3); 
% plot(omega(15,:,12,300)); hold on; 
% plot(omega_QG_RHS(15,:,12,300)); 
% legend('omega','omega_{QG}'); 
% title('Omega vs Omega_{QG} at r=0.1 and 45 deg')