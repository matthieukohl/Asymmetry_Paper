% plot lambda gcm vs. inversions

close all; clear;

list = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

strat_level = 16;

% w = 0 or w = w_gcm lower boundary condition
%low = 'w0';
low = 'w';

% take or do not take out zonal mean

averaging = 'zm';
tf_averaging = strcmp(averaging,'wzm');

lambda_turb_omega = zeros(size(r));
lambda_mode_omega = zeros(size(r));
lambda_turb_omega_QG = zeros(size(r));
lambda_turb_omega_QG_rhs = zeros(size(r));
lambda_mode_omega_QG = zeros(size(r));
lambda_mode_omega_QG_rhs = zeros(size(r));

% calculate lat variable

path_lat_turb = strcat('/disk7/mkohl/interpolation/output1/data_rescale0.1/temp_day0225h00.nc');

lat = ncread(path_lat_turb,'latitude');
lat = lat(83:111);

% mode 
tic
for tt = 1:length(r)

path = ['/net/graupel/archive1/users/mkohl/inversion/output/r_mode_simulation/rescale',list{tt},'/rms6_geostrophic_boundary_',low,'_day0560h00.mat'];

load(path,'omega');

if tf_averaging ==1

omega = omega-mean(omega,2);

end

%pholder = Lambda_weighted(omega(:,:,1:strat_level,:),lat);
pholder = Lambda_zonal(omega(:,:,1:strat_level,:));


lambda_mode_omega(tt) = mean(pholder);

clear('omega','pholder');

load(path,'omega_QG');

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

%pholder = Lambda_weighted(omega_QG(:,:,1:strat_level,:),lat);
pholder = Lambda_zonal(omega_QG(:,:,1:strat_level,:));


lambda_mode_omega_QG(tt) = mean(pholder);

clear('omega_QG','pholder');

path_rhs = ['/net/graupel/archive1/users/mkohl/inversion/output/r_mode_simulation/rescale',list{tt},'/rms6_geostrophic_boundary_',low,'_rhs_day0560h00.mat'];

load(path_rhs,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

%pholder = Lambda_weighted(omega_QG(:,:,1:strat_level,:),lat);
pholder = Lambda_zonal(omega_QG(:,:,1:strat_level,:));


lambda_mode_omega_QG_rhs(tt) = mean(pholder);

clear('omega_QG','pholder');

end
disp('mode complete')

% macroturbulent

for tt = 1:length(r)
tt

lambda_omega = 0;
lambda_omega_QG = 0;
lambda_omega_QG_rhs = 0;
    
  for jj = 1:4

path = ['/net/graupel/archive1/users/mkohl/inversion/output/r_simulation/rescale',list{tt}];

path = [path,'/rms4_geostrophic_boundary_',low,'_day',num2str(25*(8+jj)),'.mat'];

load(path,'omega');

if tf_averaging ==1

omega = omega - mean(omega,2);

end

% lambda_omega = lambda_omega + Lambda_weighted(omega(:,:,1:strat_level,:),lat);
lambda_omega = lambda_omega + Lambda_zonal(omega(:,:,1:strat_level,:));


clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

% lambda_omega_QG = lambda_omega_QG + Lambda_weighted(omega_QG(:,:,1:strat_level,:),lat);
lambda_omega_QG = lambda_omega_QG + Lambda_zonal(omega_QG(:,:,1:strat_level,:));


clear('omega_QG')

path_rhs = ['/net/graupel/archive1/users/mkohl/inversion/output/r_simulation/rescale',list{tt}];

path_rhs = [path_rhs,'/rms4_geostrophic_boundary_',low,'_rhs_day',num2str(25*(8+jj)),'.mat'];

load(path_rhs,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);

end

%lambda_omega_QG_rhs = lambda_omega_QG_rhs + Lambda_weighted(omega_QG(:,:,1:strat_level,:),lat);
lambda_omega_QG_rhs = lambda_omega_QG_rhs + Lambda_zonal(omega_QG(:,:,1:strat_level,:));


clear('omega_QG')


  end
  
lambda_turb_omega(tt) = mean(lambda_omega)/4;
lambda_turb_omega_QG(tt) = mean(lambda_omega_QG)/4;
lambda_turb_omega_QG_rhs(tt) = mean(lambda_omega_QG_rhs)/4;
end
toc

%save('figures_paper/lambda_r_w0.mat','lambda_turb_omega','lambda_turb_omega_QG','lambda_turb_omega_QG_rhs',...
 %   'lambda_mode_omega','lambda_mode_omega_QG','lambda_mode_omega_QG_rhs')

% save(['figures_paper/lambda_r_',low,'_',averaging,'.mat'],'lambda_turb_omega','lambda_turb_omega_QG','lambda_turb_omega_QG_rhs',...
%     'lambda_mode_omega','lambda_mode_omega_QG','lambda_mode_omega_QG_rhs')

% save(['figures_paper/lambda_r_',low,'_',averaging,'_rms1em5.mat'],'lambda_turb_omega','lambda_turb_omega_QG','lambda_turb_omega_QG_rhs',...
%     'lambda_mode_omega','lambda_mode_omega_QG','lambda_mode_omega_QG_rhs')

save(['figures_paper/lambda_r_',low,'_',averaging,'_zonal.mat'],'lambda_turb_omega','lambda_turb_omega_QG','lambda_turb_omega_QG_rhs',...
    'lambda_mode_omega','lambda_mode_omega_QG','lambda_mode_omega_QG_rhs')
