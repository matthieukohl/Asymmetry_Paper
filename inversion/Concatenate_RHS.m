clear;

%% RHS contribution to Skewness

path_temp = '/disk7/mkohl/interpolation/output1/data_rescale0.01/temp_day0075h00.nc';
lat_series_original = ncread(path_temp, 'latitude');
event_latspan = [79:115];
lat = lat_series_original(event_latspan);

omega_RHS_r001 = [];
omega_RHS_r01 = [];
omega_RHS_r06 = [];
omega_RHS_r1 = [];

omega = [];



for i=1:16
    
    s = load(strcat('output/rescale_0.1/sigma/equivalent_J/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    t = load(strcat('output/rescale_0.01/sigma/equivalent_J/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    u = load(strcat('output/rescale_0.6/sigma/equivalent_J/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    v = load(strcat('output/rescale_1.0/sigma/zero_BC/random/varyingf/complete/RHS_rhs_smooth0_iter30_day',num2str(25 * i),'.mat'));
    omega = cat(4,omega,s.omega);   
    omega_RHS_r001 = cat(4,omega_RHS_r001,t.omega_QG);
    omega_RHS_r01= cat(4,omega_RHS_r01,s.omega_QG);
    omega_RHS_r06 = cat(4,omega_RHS_r06,u.omega_QG);
    omega_RHS_r1 = cat(4,omega_RHS_r1,v.omega_QG);
   
end
 

omega = omega-mean(omega(:,:,:,:),2);
 
omega_RHS_r001 = omega_RHS_r001-mean(omega_RHS_r001(:,:,:,:),2);

omega_RHS_r01 = omega_RHS_r01 -mean(omega_RHS_r01(:,:,:,:),2);

omega_RHS_r06 = omega_RHS_r06-mean(omega_RHS_r06(:,:,:,:),2);

omega_RHS_r1 = omega_RHS_r1- mean(omega_RHS_r1(:,:,:,:),2);

%omega_QG_RHS = omega_QG_RHS(3:end-2,15:end-14,:,:);


              
skewness_omega = Skewness(omega);
skewness_omega_RHS_r001 = Skewness(omega_RHS_r001);
skewness_omega_RHS_r01 = Skewness(omega_RHS_r01);
skewness_omega_RHS_r06 = Skewness(omega_RHS_r06);
skewness_omega_RHS_r1 = Skewness(omega_RHS_r1);

lambda_omega = Lambda(omega);
lambda_omega_RHS_r001 = Lambda(omega_RHS_r001);
lambda_omega_RHS_r01 = Lambda(omega_RHS_r01);
lambda_omega_RHS_r06 = Lambda(omega_RHS_r06);
lambda_omega_RHS_r1 = Lambda(omega_RHS_r1);

%% Plotting


% figure(1)
% 
% plot(skewness_omega); hold on
% plot(skewness_omega_RHS_r001); hold on
% plot(skewness_omega_RHS_r01); hold on
% plot(skewness_omega_RHS_r06); hold on;
% plot(skewness_omega_RHS_r1);
% title('Skewness Vertical Velocity RHS')
% legend('omega','RHS=0.01','RHS=0.1','RHS=0.6','RHS=1.0');
% xlabel('time(quarter days)');
% ylabel('skewness')



figure(2)

t=10:1200;

plot(t,lambda_omega(10:1200)); hold on 
plot(t,lambda_omega_RHS_r001(10:1200)); hold on
plot(t,lambda_omega_RHS_r01(10:1200)); hold on
plot(t,lambda_omega_RHS_r06(10:1200)); hold on;
plot(t,lambda_omega_RHS_r1(10:1200));
title('Lambda Vertical Velocity RHS')
legend('omega','RHS=0.01','RHS=0.1','RHS=0.6','RHS=1');
xlabel('time(quarter days)');
ylabel('lambda')

saveas(gcf,'LambdaRHS.png')

