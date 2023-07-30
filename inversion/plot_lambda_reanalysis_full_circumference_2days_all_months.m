% plot lambda gcm vs. inversions reanalysis
% full circumference in the NH
% only uses 2 days per months with 6-hourly fields as a test case

close all; clear;

% reanalysis without rescale



resolution = 'coarse';
tf_resolution = strcmp(resolution,'original');

if tf_resolution==1

spec_500_omega = zeros(1,1440);
spec_500_omega_QG = zeros(1,1440);

lambda_omega_year = zeros(12,1);
lambda_omega_QG_year = zeros(12,1);

else

spec_500_omega = zeros(1,240);
spec_500_omega_QG = zeros(1,240);

lambda_omega_year = zeros(12,1);
lambda_omega_QG_year = zeros(12,1);

end

err = [];

p500 = 16;
top = 22;
counter_spectrum = 0;

for tt = 1:12

tt

path_root = '/net/graupel/archive1/users/mkohl/inversion/output/reanalysis';

if tf_resolution ==1
    
path  = [path_root,'/rms1em4_geostrophic_boundary_w_without_rescale_whole_circumference_t',list{tt},'.mat'];

else
    
path  = [path_root,'/rms1em4_boundary_w_without_rescale_whole_circumference_geostrophic_res_6_2days_for_all_months.mat'];

end

load(path,'lat'); load(path,'lon')

% calculate Lambda

load(path,'omega');
load(path,'omega_QG');

[a,lonm180] = min(abs(lon+180));
[a,lonm120] = min(abs(lon+120));

[a,lonm60] = min(abs(lon+60));
[a,lonm10] = min(abs(lon+10));

[a,lon130] = min(abs(lon-130));
[a,lon180] = min(abs(lon-180));

%ind_lon = [lonm180:lonm120,lonm60:lonm10,lon130:lon180];
ind_lon = [1:length(lon)];

counter = 0;
lambda_omega = 0;
lambda_omega_QG = 0;

for t = 1:8
    
ind = (tt-1)*8 + t;

lambda_omega = lambda_omega + Lambda_weighted(omega(:,ind_lon,p500,ind),lat);

lambda_omega_QG = lambda_omega_QG + Lambda_weighted(omega_QG(:,ind_lon,p500,ind),lat);

counter = counter + 1;

end

% calculate the spectra

for t = 1:8
 
ind = (tt-1)*8 + t;
    
[k_500_new,spec_500_true_new] = Spectrum2(omega(:,:,p500,ind),lat,lon);
[k_500_new,spec_500_inverted_new] = Spectrum2(omega_QG(:,:,p500,ind),lat,lon);
spec_500_omega = spec_500_omega + spec_500_true_new;
spec_500_omega_QG = spec_500_omega_QG + spec_500_inverted_new;

counter_spectrum = counter_spectrum + 1;
% % inversions
% 
% R = 6371000.0; 
% R_eff = R*cosd(45);
% N = length(lon);
% k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
% Ld = 1e6/(2*sqrt(2));
% k = k*Ld;
% 
% figure
% loglog(k(1:N/2),spec_500_true_new(1:N/2)/max(spec_500_true_new(1:N/2))); hold on;
% loglog(k(1:N/2),spec_500_inverted_new(1:N/2)/max(spec_500_inverted_new(1:N/2))); hold on;

end

lambda_omega_year(tt) = lambda_omega/counter;
lambda_omega_QG_year(tt) = lambda_omega_QG/counter;

load(path,'err_rms');

ind1 = (tt-1)*8 + 1;
ind2 = (tt-1)*8 + 8;

err = cat(1,err,err_rms(ind1:ind2));


end

% yearly avaeraged spectra

spec_500_omega = spec_500_omega/counter_spectrum;
spec_500_omega_QG = spec_500_omega_QG/counter_spectrum;


% plot lambda seasonal cycle vor omega vs omega QG

figure
plot(lambda_omega_year,'b','Linewidth',1.4); hold on;
plot(lambda_omega_QG_year,'b--','Linewidth',1.4);
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Asymmetry Parameter $\lambda$');
set(gca,'Fontsize',12)
legend('$\omega$','$\omega$ QG','Location','NorthWest'); legend boxoff;


% plot yearly average lambda for omega vs omega QG

lambda = [mean(lambda_omega_year),mean(lambda_omega_QG_year)];

figure

labels={'omega','QG'};

bar(lambda)
set(gca,'xticklabel',labels);
ylabel('$\lambda$')
title('\rm Reanalysis Comparison')
ylim([0 0.9])
set(gca,'FontSize',12)

% plot spectra

R = 6371000.0; 
R_eff = R*cosd(45);
N = length(lon);
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
%Ld = 1e6/(2*sqrt(2)); 1000km estimate
Ld = 4.67*1e5;
k = k*Ld;

figure
loglog(k(1:N/2),spec_500_omega(1:N/2)/max(spec_500_omega(1:N/2))); hold on;
loglog(k(1:N/2),spec_500_omega_QG(1:N/2)/max(spec_500_omega_QG(1:N/2))); hold on;
loglog(k(1:N/2),1e1*k(1:N/2).^(-1),'k--'); hold on;
xlabel('Wavenumber k');
ylabel('Spec')
legend('true','inverted','k^{-1}','Location','SouthWest'); legend boxoff
title('\rm North Atlantic (500hPa)')

k_centroid_spec_omega = sum(k(1:N/2).*spec_500_omega(1:N/2)/max(spec_500_omega(1:N/2)))/sum(spec_500_omega(1:N/2)/max(spec_500_omega(1:N/2)));
k_centroid_spec_omega_QG = sum(k(1:N/2).*spec_500_omega_QG(1:N/2)/max(spec_500_omega_QG(1:N/2)))/sum(spec_500_omega_QG(1:N/2)/max(spec_500_omega_QG(1:N/2)));


k_centroid = [k_centroid_spec_omega,k_centroid_spec_omega_QG];

figure

labels={'omega','QG'};

bar(k_centroid)
set(gca,'xticklabel',labels);
ylabel('k')
title('\rm Wavenumber')
set(gca,'FontSize',12)

% wavenumber cut-off at L = 500km

k500 = 2*pi/(500*1e3)*Ld;

[a,ind500] = min(abs(k(1:N/2)-k500));

k_centroid_spec_omega = sum(k(1:ind500).*spec_500_omega(1:ind500)/max(spec_500_omega(1:ind500)))/sum(spec_500_omega(1:ind500)/max(spec_500_omega(1:ind500)));
k_centroid_spec_omega_QG = sum(k(1:ind500).*spec_500_omega_QG(1:ind500)/max(spec_500_omega_QG(1:ind500)))/sum(spec_500_omega_QG(1:ind500)/max(spec_500_omega_QG(1:ind500)));


k_centroid = [k_centroid_spec_omega,k_centroid_spec_omega_QG];

figure

labels={'omega','QG'};

bar(k_centroid)
set(gca,'xticklabel',labels);
ylabel('k')
title('\rm Wavenumber (cut-off at 500km)')
set(gca,'FontSize',12)



% plot one of the inversions

cmax = 1;
nint = 20;
cint = cmax/nint;

figure('Renderer', 'painters', 'Position', [10 10 1200 600])

subplot(2,1,1)

contourf(-omega(:,:,p500,end),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega')

subplot(2,1,2)

contourf(-omega_QG(:,:,p500,end),-cmax:cint:cmax,'EdgeColor','none'); hold on;
colormap(redblue(nint));
colorbar;
caxis([-cmax cmax])
title('\rm Omega Inverted')

% plot the errors

figure

plot(err,'b'); hold on;
plot(1e-4*ones(size(err)),'k');
xlabel('inversions')
ylabel('rms err')
title('\rm without rescale')