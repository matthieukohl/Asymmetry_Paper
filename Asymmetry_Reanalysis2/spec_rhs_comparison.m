% compare the rhs-spectra in gcm, qg, and era5

close all; 
clear;

% load gcm simulation

load('/disk7/mkohl/inversion/spec_rhs_gcm_r001.mat');

k_gcm = k*1e6/(2*sqrt(2));
spec_rhs_gcm = spec_500_r(1,:); % select the r=0.01 run
spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);

clear('k','spec_500_r');

% load qg simulation

load('spec_rhs_qg_r001.mat')

k_qg = kx;
spec_rhs_qg = power_rhs_kx;

spec_rhs_qg = spec_rhs_qg/max(spec_rhs_qg);

[a,ind_qg] = min(abs(kx-36));

clear('kx','power_rhs_kx');

% load era5 

load('spec_rhs_era5.mat');

k_era5 = k;
spec_rhs_era5 = spec_rhs_winter;

spec_rhs_era5 = spec_rhs_era5/max(spec_rhs_era5);

[a,ind_era5] = min(abs(k-38));

clear('k','spec_rhs_winter');

% make the plot

figure(1)
loglog(k_gcm,spec_rhs_gcm,'b'); hold on;
loglog(k_qg(1:ind_qg),spec_rhs_qg(1:ind_qg),'r'); hold on;
loglog(k_era5(1:ind_era5),spec_rhs_era5(1:ind_era5),'k');
xlabel('Wavenumber k')
ylabel('Amplitude (normalized)')
legend('GCM','QG','ERA5'); legend boxoff
title('\rm 1-D RHS Spectrum');
set(gca,'FontSize',12);
ylim([1e-3 1e1])


% fit slopes

figure(2)
loglog(k_gcm,spec_rhs_gcm,'b'); hold on;
ylim([1e-3 1e1])

[t,ind] = max(abs(squeeze(spec_rhs_gcm)));

k1 = k_gcm(2);
k2 = k_gcm(ind);

% calculate power spectrum of slope alpha over k span k1 to k2

Lx = 2*pi/k1;
N = k2*Lx/pi;

xx = linspace(0,Lx,N);

k = 2*pi/(Lx)*[0:N/2,-N/2+1:-1];

loglog(k,k.^1/100,'r')



