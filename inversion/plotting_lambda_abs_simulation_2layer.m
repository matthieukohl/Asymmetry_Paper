% plot lambda gcm vs. inversions abs simulations
close all; clear;

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

%strat_level = [8,9,10,10,11,11,12,12,13,13,14,14,15,16,18];
%strat_level = strat_level + 5;

load('figures_paper/trop_height.mat');
strat_level = trop_height;

% boundary condition w = 0 or w = w_gcm

%low = 'w0';
low = 'w';

% take or do not take out zonal mean

averaging = 'zm';
tf_averaging = strcmp(averaging,'wzm');

lambda_turb_omega = zeros(length(list),1);
lambda_mode_omega = zeros(length(list),1);
lambda_turb_omega_QG = zeros(length(list),1);
lambda_turb_omega_QG_rhs = zeros(length(list),1);
lambda_mode_omega_QG = zeros(length(list),1);
lambda_mode_omega_QG_rhs = zeros(length(list),1);

% calculate lat variable

path_lat_turb = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.1/temp_day0225h00.nc');

lat = ncread(path_lat_turb,'latitude');
lat = lat(83:111);
lat = lat(11:26);

% mode 
tic
for tt = 1:length(list)

path = ['/net/graupel/archive1/users/mkohl/inversion/output/abs_mode_simulation/',list{tt},'/rms15_geostrophic_',low,'_day0560h00_two_layer_alternative.mat'];

load(path,'omega');

if tf_averaging==1

omega = omega-mean(omega,2);

end

pholder = Lambda_weighted(omega(11:26,:,2,:),lat);
%pholder = Lambda(omega(11:26,:,5:strat_level(tt),:));

lambda_mode_omega(tt) = mean(pholder);

clear('omega','pholder');

load(path,'omega_QG');

if tf_averaging==1

omega_QG = omega_QG - mean(omega_QG,2);

end

pholder = Lambda_weighted(omega_QG(11:26,:,2,:),lat);
%pholder = Lambda(omega_QG(11:26,:,5:strat_level(tt),:));


lambda_mode_omega_QG(tt) = mean(pholder);

clear('omega_QG','pholder');

% path_rhs = ['/net/graupel/archive1/users/mkohl/inversion/output/abs_mode_simulation/',list{tt},'/rms15_geostrophic_',low,'_rhs_day0560h00.mat'];
% 
% load(path_rhs,'omega_QG');
% 
% if tf_averaging ==1
% 
% omega_QG = omega_QG - mean(omega_QG,2);
% 
% end
% 
% pholder = Lambda_weighted(omega_QG(11:26,:,1:strat_level(tt),:),lat);
% %pholder = Lambda(omega_QG(11:26,:,5:strat_level(tt),:));
% 
% 
% lambda_mode_omega_QG_rhs(tt) = mean(pholder);
% 
% clear('omega_QG','pholder');

end
disp('mode complete')

% macroturbulent

for tt = 1:length(list)
tt

lambda_omega = 0;
lambda_omega_QG = 0;
lambda_omega_QG_rhs = 0;

counter = 0;
    
  for jj = 1:6

path = ['/net/graupel/archive1/users/mkohl/inversion/output/abs_simulation/',list{tt}];

path = [path,'/rms1em4_geostrophic_boundary_',low,'_day_',num2str(5*(28+jj)),'_2layer_alternative.mat'];

load(path,'omega');

if tf_averaging ==1

omega = omega - mean(omega,2);

end

lambda_omega = lambda_omega + Lambda_weighted(omega(11:26,:,2,:),lat);

clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

lambda_omega_QG = lambda_omega_QG + Lambda_weighted(omega_QG(11:26,:,2,:),lat);

clear('omega_QG')

% path_rhs = ['/net/graupel/archive1/users/mkohl/inversion/output/abs_simulation/',list{tt}];
% 
% path_rhs = [path_rhs,'/rms1em4_geostrophic_boundary_',low,'_rhs_day',num2str(5*(28+jj)),'.mat'];
% 
% load(path_rhs,'omega_QG');
% 
% if tf_averaging ==1
%     
% omega_QG = omega_QG - mean(omega_QG,2);
% 
% end
% 
% lambda_omega_QG_rhs = lambda_omega_QG_rhs + Lambda_weighted(omega_QG(11:26,:,1:strat_level(tt),:),lat);

counter = counter + 1;

  end
  
lambda_turb_omega(tt) = mean(lambda_omega)/counter;
lambda_turb_omega_QG(tt) = mean(lambda_omega_QG)/counter;
%lambda_turb_omega_QG_rhs(tt) = mean(lambda_omega_QG_rhs)/jj;
  
end
toc

% save(['figures_paper/test_lambda_abs_',low,'_',averaging,'.mat'],'lambda_turb_omega','lambda_turb_omega_QG','lambda_turb_omega_QG_rhs',...
%     'lambda_mode_omega','lambda_mode_omega_QG','lambda_mode_omega_QG_rhs')

% save(['figures_paper/lambda_abs_',low,'_',averaging,'_rms1em4_2layer_alternative.mat'],'lambda_turb_omega','lambda_turb_omega_QG',...
%     'lambda_mode_omega','lambda_mode_omega_QG')

% % boundary condition w = 0 or w = w_gcm
% 
% %low = 'w0';
% low = 'w';
% 
% % take or do not take out zonal mean
% 
% averaging = 'zm';
% tf_averaging = strcmp(averaging,'wzm');
% 
% load(['figures_paper/lambda_abs_',low,'_',averaging,'_rms1em4_2layer.mat'],'lambda_turb_omega','lambda_turb_omega_QG',...
%     'lambda_mode_omega','lambda_mode_omega_QG')


surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

figure(1)
%x0=10; y0=10; width=1000; height=400; set(gcf,'position',[x0,y0,width,height])

plot(surf_temp,lambda_turb_omega,'b','linewidth',1.2); hold on;
plot(surf_temp,lambda_turb_omega_QG,'b--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega,'r','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega_QG,'r--');
xlabel('Global Surface Temperature (K)');
ylabel('$\lambda$');
set(gca,'FontSize',12);
title('\rm 2-Layer Inversions') 
ylim([0.5 1])
set(gca,'box','off')

%saveas(gcf,['figures_paper/inversions_abs_',low],'epsc')
