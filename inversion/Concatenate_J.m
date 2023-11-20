clear;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0075h00.nc';
lat_series_original = ncread(path_temp, 'latitude');
event_latspan = [79:115];
lat = lat_series_original(event_latspan);

omega_QG_RHS = [];
omega_QG_LHS = [];
omega = [];

for i=1:16
    
    s = load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/RHS_rhs_smooth0_iter10_day',num2str(25 * i),'.mat'));
    %t = load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/LHS_rhs_smooth0_iter10_day',num2str(25 * i),'.mat'));
    omega = cat(4,omega,s.omega);  
    omega_QG_RHS=cat(4,omega_QG_RHS,s.omega_QG);
    %omega_QG_LHS = cat(4,omega_QG_LHS,t.omega_QG);
end
 

omega = omega-mean(omega(:,:,:,:),2);
 
omega_QG_RHS = omega_QG_RHS-mean(omega_QG_RHS(:,:,:,:),2);

omega_QG_LHS = omega_QG_LHS -mean(omega_QG_LHS(:,:,:,:),2);

%omega_QG_RHS = omega_QG_RHS(3:end-2,15:end-14,:,:);


              
skewness_omega = Skewness(omega);
skewness_omegaQG_RHS = Skewness(omega_QG_RHS);
skewness_omega_QG_LHS = Skewness(omega_QG_LHS);

lambda_omega = Lambda(omega);
lambda_omegaQG_RHS = Lambda(omega_QG_RHS);
lambda_omegaQG_LHS = Lambda(omega_QG_LHS);


%% Total Contribution to Skewness

omega_QG_total= [];
sigma_new= [];
sigma_new_smoothed= [];

for i=1:16
    
    t= load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/random/Total_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    omega_QG_total=cat(4,omega_QG_total,t.omega_QG);
        
end


omega = omega-mean(omega(:,:,:,:),2);

omega_QG_total = omega_QG_total-mean(omega_QG_total(:,:,:,:),2);

%omega_QG_total = omega_QG_total(3:end-2,15:end-14,:,:);


skewness_omegaQG_total = Skewness(omega_QG_total);
lambda_omegaQG_total = Lambda(omega_QG_total);




%% Plotting


figure(3)

plot(skewness_omega); hold on
plot(skewness_omegaQG_RHS); hold on
%plot(skewness_omega_QG_LHS); hold on
plot(skewness_omegaQG_total);
title('Skewness Vertical Velocity (r=0.1)')
legend('omega','omegaQG_{RHS}','omega_QG_{LHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('skewness')



figure(4)

plot(lambda_omega); hold on 
plot(lambda_omegaQG_RHS); hold on
%plot(lambda_omegaQG_LHS); hold on
plot(lambda_omegaQG_total);
title('Lambda Vertical Velocity (r=0.1)')
legend('omega','omegaQG_{RHS}','omegaQG_{LHS}','omegaQG_{total}');
xlabel('time(quarter days)');
ylabel('lambda')



% figure(3); 
% plot(omega(15,:,12,200)); hold on; 
% plot(omega_QG_total(15,:,12,200)); 
% legend('omega','omega_{QG}'); 
% title('Omega vs Omega_{QG} at r=0.1 and 45 deg')
% 
% 
% 
% %% Averages
% 
% m = squeeze(mean(mean(omega_QG_total(:,:,:,:),1),2));
% n = squeeze(mean(mean(omega_QG_RHS(:,:,:,:),1),2));
% o = squeeze(mean(mean(omega(:,:,:,:),1),2));
% 
% 
% m = mean(m(:,:),1);
% n = mean(n(:,:),1);
% o = mean(o(:,:),1);
% 
% s = weighted_avg(omega_QG_total,lat);
% t = weighted_avg(omega_QG_RHS,lat);
% u = weighted_avg(omega,lat);
% 
% figure(4); 
% 
% plot(m); hold on; 
% plot(n); hold on;
% plot(o)
% 
% title('Averages')
% legend('QG_{total}','QG_{RHS}','Omega')
% 
% figure(5)
% 
% plot(s); hold on; 
% plot(t); hold on;
% plot(u);
% title('Weighted Averages')
% legend('QG_{total}','QG_{RHS}','Omega')



% figure(3)
% 
% for i=1:length(omega,4)
% contour(omega_QG_LHS(:,:,12,i)); 
% pause(0.2)
% end

% figure(1)
% 
% histogram(omega_QG_total(:,:,12,:),200); hold on; 
% title('PDF-Omega_{Total}')
% 
% figure(2); 
% histogram(omega_QG_RHS(:,:,12,:),200); hold on;
% title('PDF-Omega_{RHS}')
% 
% figure(3);
% histogram(omega(:,:,12,:),200); hold on;
% title('PDF-Omega')



