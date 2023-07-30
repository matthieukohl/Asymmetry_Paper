% plot lambda gcm vs. inversions abs simulations

close all; clear;

% boundary condition w = 0 or w = w_gcm

%low = 'w0';
low = 'w';

% take or do not take out zonal mean

averaging = 'wzm';
tf_averaging = strcmp(averaging,'wzm');

% load levels

path_level = '/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_atlantic_patch_inversion_with_uv.nc';

level = ncread(path_level,'level');
level = double(level)*100; level = flip(level);
top = 22;
[a,p500] = min(abs(level-500*1e2));


% reanalysis without rescale

% high resolution
path = '/net/graupel/archive1/users/mkohl/inversion/output/reanalysis/rms1em4_geostrophic_boundary_w_without_rescale.mat';

load(path,'lat');

load(path,'omega');


if tf_averaging ==1

omega = omega - mean(omega,2);
%omega = omega - mean(mean(omega,1),2);

end

lambda_omega = Lambda_weighted(omega(:,:,1:top,:),lat);
lambda_omega = mean(lambda_omega);

% lambda_omega = zeros(120,1);
% 
% for tt = 1:120
% 
% lambda_omega(tt) = Lambda_weighted(omega(:,:,p500,tt),lat); 
% end
% lambda_omega = mean(lambda_omega);

clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);
%omega_QG = omega_QG - mean(mean(omega_QG,1),2);
end

lambda_omega_QG = Lambda_weighted(omega_QG(:,:,1:top,:),lat);
lambda_omega_QG = mean(lambda_omega_QG);

% lambda_omega_QG = zeros(120,1);
% for tt = 1:120
% 
% lambda_omega_QG(tt) = Lambda_weighted(omega_QG(:,:,p500,tt),lat); 
% end
% lambda_omega_QG = mean(lambda_omega_QG);



clear('omega_QG')

load(path,'err_rms')
err_true = err_rms;

clear('err_rms')

% coarse res

path = '/net/graupel/archive1/users/mkohl/inversion/output/reanalysis/rms1em4_geostrophic_boundary_w_without_rescale_true_coarse_resolution.mat';

load(path,'omega');
load(path,'lat');


if tf_averaging ==1

omega = omega - mean(omega,2);
%omega = omega - mean(mean(omega,1),2);
end

lambda_omega_coarse = Lambda_weighted(omega(:,:,1:top,:),lat);
lambda_omega_coarse = mean(lambda_omega_coarse);

% lambda_omega_coarse = zeros(120,1);
% for tt = 1:120
% 
% lambda_omega_coarse(tt) = Lambda_weighted(omega(:,:,p500,tt),lat); 
% end
% lambda_omega_coarse = mean(lambda_omega_coarse);
% 

clear('omega')

load(path,'omega_QG')

if tf_averaging ==1

omega_QG = omega_QG - mean(omega_QG,2);
%omega_QG = omega_QG - mean(mean(omega_QG,1),2);
end

lambda_omega_QG_coarse = Lambda_weighted(omega_QG(:,:,1:top,:),lat);
lambda_omega_QG_coarse = mean(lambda_omega_QG_coarse);

% lambda_omega_QG_coarse = zeros(120,1);
% for tt = 1:120
% 
% lambda_omega_QG_coarse(tt) = Lambda_weighted(omega_QG(:,:,p500,tt),lat); 
% end
% lambda_omega_QG_coarse = mean(lambda_omega_QG_coarse);


clear('omega_QG')

load(path,'err_rms')
err_coarse = err_rms;

clear('err_rms')

%
lambda = [lambda_omega,lambda_omega_QG,lambda_omega_coarse,lambda_omega_QG_coarse];

figure

labels={'omega','QG','omega (coarse)','QG (coarse)'};

bar(lambda)
set(gca,'xticklabel',labels);
ylabel('$\lambda$')
title('\rm Reanalysis Comparison Coarse-Grained')
ylim([0 0.9])

% plot the errors

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)
plot(err_true,'b'); hold on;
plot(1e-4*ones(size(err_true)),'k');
xlabel('inversions')
ylabel('rms err')
title('\rm without rescale')

subplot(1,2,2)
plot(err_coarse,'r'); hold on;
plot(1e-4*ones(size(err_coarse)),'k');
xlabel('inversions')
ylabel('rms err')
title('\rm with rescale')

load(path,'omega','omega_QG','lat')

% plot
cmax = 1;
nint = 20;
cint = cmax/nint;
[a,p500] = min(abs(level-50000));

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

subplot(2,1,1)

contourf(-omega(:,:,p500,100),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega')

subplot(2,1,2)

contourf(-omega_QG(:,:,p500,100),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega Inverted')

% % calculate the spectrum of omega vs omega QG
% 
% R = 6371000.0; % Planetary Radius
% 
% [a,p500] = min(abs(level-500*1e2));
% 
% spec_500_true = zeros(1,161);
% spec_500_inverted = zeros(1,161);
% counter = 0;
% 
% for tt = 1:120
% [k_500_new,spec_500_true_new] = Spectrum2(omega(:,:,p500,tt),lat,lon);
% [k_500_new,spec_500_inverted_new] = Spectrum2(omega_QG(:,:,p500,tt),lat,lon);
% spec_500_true = spec_500_true + spec_500_true_new;
% spec_500_inverted = spec_500_inverted + spec_500_inverted_new;
% counter = counter + 1;
% end
% 
% spec_500_true = spec_500_true/counter; 
% spec_500_inverted = spec_500_inverted/counter;
% 
% R_eff = R*cosd(45);
% %N = length(lon);
% N = size(omega,2);
% k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
% Ld = 1e6/(2*sqrt(2));
% k = k*Ld;
% 
% figure
% loglog(k(1:N/2),spec_500_true(1:N/2)/max(spec_500_true(1:N/2))); hold on;
% loglog(k(1:N/2),spec_500_inverted(1:N/2)/max(spec_500_inverted(1:N/2))); hold on;
% loglog(k(1:N/2),1e1*k(1:N/2).^(-1),'k--'); hold on;
% xlabel('Wavenumber k');
% ylabel('Spec')
% legend('true','inverted','k^{-1}','Location','SouthWest'); legend boxoff
% title('\rm North Atlantic (500hPa)')

% test