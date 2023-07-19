clear;

%% Concatenation Variables

omega_QG = [];
omega = [];


% for i=1:12
%     
%     s= load(strcat('output/rescale_0.6/r_0.6_Totalcontribution_day',num2str(25 * i),'.mat'));
%     
%     omega = cat(4,omega,s.omega);
%     omega_QG = cat(4,omega_QG,s.omega_QG);
% end

for i=1:12
    
    s= load(strcat('output/rescale_0.6/sigma/equivalent_J/zero_BC/Total_rhs_smooth0_iter10_day',num2str(25 * i),'.mat'));
    
    omega = cat(4,omega,s.omega);
    omega_QG = cat(4,omega_QG,s.omega_QG);
end


omega=omega-mean(omega(:,:,:,:),2);

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
% title('Skewness Vertical Velocity (r=0.6)')
% legend('omega','omega_{QG}');
% xlabel('time(quarter days)');
% ylabel('skewness')
% 
% savefig('output/rescale_0.6/skewness.fig');


figure(2)

t=200:1200;

plot(t,lambda_omega(200:end)); hold on 
plot(t,lambda_omegaQG(200:end)); hold on
title('Lambda Vertical Velocity (r=0.6)')
legend('omega','omega_{QG}');
xlabel('time(quarter days)');
ylabel('lambda')

saveas(gcf,'Lambda06.png')

%% Averages

% m = squeeze(mean(mean(omega_QG(:,:,:,:),1),2));
% n = squeeze(mean(mean(omega(:,:,:,:),1),2));
% 
% 
% m = mean(m(:,:),1);
% n = mean(n(:,:),1);
% 
% % s = weighted_avg(omega_QG,lat);
% % t = weighted_avg(omega,lat);
% 
% figure(4); 
% 
% plot(m); hold on; 
% plot(n);
% title('Averages')
% legend('QG_{total}','Omega')

% figure(5)
% 
% plot(s); hold on; 
% plot(t);
% title('Weighted Averages')
% legend('QG_{total}','Omega}')



