clear;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0075h00.nc';
lat_series_original = ncread(path_temp, 'latitude');
event_latspan = [83:111];
lat = lat_series_original(event_latspan);

omega_QG_RHS = [];
omega = [];

for i=1:3
    
    s = load(strcat('output/rescale_0.1/sigma/relax_mean/RHS_rhs_smooth0_sigma_iter10_day',num2str(25 * i),'.mat'));
    omega = cat(4,omega,s.omega);  
    omega_QG_RHS=cat(4,omega_QG_RHS,s.omega_QG);
end
    

omega = omega-mean(omega(:,:,:,:),2);
 
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);

skewness_omega = Skewness(omega);
skewness_omegaQG_RHS = Skewness(omega_QG_RHS);

lambda_omega = Lambda(omega);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS);


%% Total Contribution to Skewness

omega_QG_total= [];
sigma_new= [];
sigma_new_smoothed= [];

for i=1:3
    
    t= load(strcat('output/rescale_0.1/sigma/relax_mean/smooth_cross/Total_rhs_smooth0_sigma_smooth4.0_iter10_day',num2str(25 * i),'.mat'));
    
    omega_QG_total=cat(4,omega_QG_total,t.omega_QG);
    sigma_new=cat(4,sigma_new,t.sigma_new);
    sigma_new_smoothed=cat(4,sigma_new_smoothed,t.sigma_new_smoothed);
    
end

omega_QG_total = omega_QG_total-mean(omega_QG_total(:,:,:,:),2);

skewness_omegaQG_total = Skewness(omega_QG_total);
lambda_omegaQG_total = Lambda(omega_QG_total);




%% Plotting


figure(1)

plot(skewness_omega); hold on
plot(skewness_omegaQG_RHS); hold on
plot(skewness_omegaQG_total);
title('Skewness Vertical Velocity (r=0.1)')
legend('omega','omegaQG_{RHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('skewness')
axis([0 300 0 4])
savefig('output/rescale_0.1/skewness.fig');

figure(2)

plot(lambda_omega); hold on 
plot(lambda_omegaQG_RHS); hold on
plot(lambda_omegaQG_total);
title('Lambda Vertical Velocity (r=0.1)')
legend('omega','omegaQG_{RHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('lambda')

savefig('output/rescale_0.1/lambda.fig');

figure(3); 
plot(omega(15,:,12,300)); hold on; 
plot(omega_QG_total(15,:,12,300)); 
legend('omega','omega_{QG}'); 
title('Omega vs Omega_{QG} at r=0.1 and 45 deg')



%% Averages

m = squeeze(mean(mean(omega_QG_total(:,:,:,:),1),2));
n = squeeze(mean(mean(omega_QG_RHS(:,:,:,:),1),2));


m = mean(m(:,:),1);
n = mean(n(:,:),1);

s = weighted_avg(omega_QG_total,lat);
t = weighted_avg(omega_QG_RHS,lat);

figure(4); 

plot(m); hold on; 
plot(n);
title('Averages')
legend('QG_{total}','QG_{RHS}')

figure(5)

plot(s); hold on; 
plot(t);
title('Weighted Averages')
legend('QG_{total}','QG_{RHS}')

figure(6)
plot(sigma_new(12,:,12,300)); hold on;
plot(sigma_new_smoothed(12,:,12,300));
title('sigma')
legend('full','smoothed')

% m=mean(mean(omega_RHS_prime.^3,1),2)./(mean(mean(omega_RHS_prime.^2),1),2).^(3/2);
% 
% av = mean(omega_QG_RHS(15,:,12,300));
% 
% prime = omega_QG_RHS(15,:,12,300)-av;
% 
% skew = mean(prime.^3)./(mean(prime.^2)).^(3/2);






