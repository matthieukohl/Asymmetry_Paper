% plot lambda gcm vs. inversions abs simulations
close all; clear;

list = {'abs4.0_T85','abs4.0_T170'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

%strat_level = [8,9,10,10,11,11,12,12,13,13,14,14,15,16,18];
%strat_level = strat_level + 5;

load('figures_paper/trop_height.mat');
strat_level = trop_height(14); % abs4.0


% boundary condition w = 0 or w = w_gcm

%low = 'w0';
low = 'w';

% take or do not take out zonal mean

averaging = 'wzm';
tf_averaging = strcmp(averaging,'wzm');

lambda_turb_omega = zeros(length(list),1);
lambda_turb_omega_QG = zeros(length(list),1);

% calculate lat variable

path_lat_turb = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.1/temp_day0225h00.nc');

lat = ncread(path_lat_turb,'latitude');
lat = lat(83:111);
lat = lat(11:26);

% macroturbulent T85

err_T85 = [];

tic
for tt = 1:1
tt

lambda_omega = 0;
lambda_omega_QG = 0;

    
for jj = 1:6

path = ['/net/graupel/archive1/users/mkohl/inversion/output/abs_simulation/',list{tt}];

path = [path,'/rms1em4_geostrophic_boundary_',low,'_day',num2str(5*(28+jj)),'.mat'];

load(path,'omega');

if tf_averaging ==1

omega = omega - mean(omega,2);

end

lambda_omega = lambda_omega + Lambda_weighted(omega(11:26,:,1:strat_level,:),lat);

clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

lambda_omega_QG = lambda_omega_QG + Lambda_weighted(omega_QG(11:26,:,1:strat_level,:),lat);

clear('omega_QG')

load(path,'err_rms');

err_T85 = cat(1,err_T85,err_rms);

end
  
lambda_turb_omega(tt) = mean(lambda_omega)/jj;
lambda_turb_omega_QG(tt) = mean(lambda_omega_QG)/jj;

  
end


% macroturbulent T170
i = 22; tt = 2;

path_lat_turb = strcat('/net/graupel/archive1/users/mkohl/interpolation/',list{tt},'/temp_day',num2str(5 *(521+i)),'h00.nc');
lat = ncread(path_lat_turb,'latitude');
lat = lat(186:215);

err_T170 = [];

for tt = 2:2
tt

lambda_omega = 0;
lambda_omega_QG = 0;

counter = 0;
err = 0;
for i= 22:59
    i

path = strcat('/net/graupel/archive1/users/mkohl/inversion/output/abs_simulation/',list{tt},'/rms1em4_geostrophic_boundary_w_day',num2str(5 *(521+i)),'_T170.mat');

load(path,'omega');

if tf_averaging ==1

omega = omega - mean(omega,2);

end

lambda_omega = lambda_omega + Lambda_weighted(omega(:,:,1:strat_level,:),lat);

clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

lambda_omega_QG = lambda_omega_QG + Lambda_weighted(omega_QG(:,:,1:strat_level,:),lat);

clear('omega_QG')

load(path,'err_rms');

err_T170 = cat(1,err_T170,err_rms);

counter = counter + 1;

end
  
lambda_turb_omega(tt) = mean(lambda_omega)/counter;
lambda_turb_omega_QG(tt) = mean(lambda_omega_QG)/counter;

  
end
toc

% compare lambda

figure(1)

labels = {'GCM T85','QG T85','GCM T170','QG T170','diff T85','diff T170'};

lambda_turb = [lambda_turb_omega(1),lambda_turb_omega_QG(1),...
    lambda_turb_omega(2),lambda_turb_omega_QG(2),...
    lambda_turb_omega_QG(1)-lambda_turb_omega(1),lambda_turb_omega_QG(2)-lambda_turb_omega(2)];

bar(lambda_turb)
set(gca,'xticklabel',labels');
ylabel('\lambda')
title('\rm abs 4.0')

% plot the errors

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)
plot(err_T85,'b'); hold on;
plot(1e-4*ones(size(err_T85)),'k');
xlabel('inversions')
ylabel('rms err')
title('\rm T85')

subplot(1,2,2)
plot(err_T170,'r'); hold on;
plot(1e-4*ones(size(err_T170)),'k');
xlabel('inversions')
ylabel('rms err')
title('\rm T170')






