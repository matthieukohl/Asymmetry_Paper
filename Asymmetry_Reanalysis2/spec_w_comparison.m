% compare the w spectrum in gcm, qg and era5

close all; clear;

% load gcm simulation

load('/disk7/mkohl/inversion/spec_w_gcm_r001.mat');

R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

k_gcm = k*1e6/(2*sqrt(2));
%k_gcm = k*5.79*1e5; based on real N2
spec_w_gcm = spec_500_r(1,:); % select the r=0.01 run
spec_w_gcm = spec_w_gcm/max(spec_w_gcm);

spec_w_gcm_r01 = spec_500_r(4,:);
spec_w_gcm_r01 = spec_w_gcm_r01/max(spec_w_gcm_r01);

spec_w_gcm_r04 = spec_500_r(7,:);
spec_w_gcm_r04 = spec_w_gcm_r04/max(spec_w_gcm_r04);

clear('k','spec_500_r');

% load the qg simulation r=0.01

load('spec_w_qg_r001.mat')

k_qg = kx;
spec_w_qg = power_rhs_kx;

spec_w_qg = spec_w_qg/max(spec_w_qg);

[a,ind_qg] = min(abs(kx-36));

clear('kx','power_rhs_kx');

% load the qg simulation r=0.4

load('spec_w_qg_r04.mat')

k_qg_r04 = kx;

spec_w_qg_r04 = power_rhs_kx;

spec_w_qg_r04 = spec_w_qg_r04/max(spec_w_qg_r04);

%[a,ind_qg] = min(abs(kx-36));

clear('kx','power_rhs_kx');

% load the qg simulation r=1

load('spec_w_qg_r1.mat')

k_qg_r1 = kx;

spec_w_qg_r1 = power_rhs_kx;

spec_w_qg_r1 = spec_w_qg_r1/max(spec_w_qg_r1);

%[a,ind_qg] = min(abs(kx-36));

clear('kx','power_rhs_kx');

% load the qg simulation 128/128 resolution

load('spec_w_qg_r001_128_res.mat')

k_qg_low_res = kx;
spec_w_qg_low_res = power_rhs_kx;

spec_w_qg_low_res = spec_w_qg_low_res/max(spec_w_qg_low_res);

%[a,ind_qg] = min(abs(kx-36));

clear('kx','power_rhs_kx');

% load the era5 simulation

load('data_theory/theory_NH_oct13th.mat'); 


load('lon.mat');

R = 6371000.0; % Planetary Radius 
R_eff = cosd(45)*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);

k_era5 = k*1e6/(2*sqrt(2));


[a,ind] = min(abs(k));

spec_w_era5 = mean(spec_NH_500,2);
spec_w_era5 = spec_w_era5/max(spec_w_era5);

clear('k','spec_NH_500');

% make the plot

[a,indLL] = min(abs(k_gcm-0.5));
[a,indL] = min(abs(k_gcm-1.4));

figure(1)
loglog(k_gcm,spec_w_gcm,'b','Linewidth',1.4); hold on;
%loglog(k_gcm,spec_w_gcm_r01,'b:','Linewidth',1.4); hold on;
loglog(k_qg(1:ind_qg),spec_w_qg(1:ind_qg),'r','Linewidth',1.4); hold on;
%loglog(k_qg(1:ind_qg),spec_w_qg_r04(1:ind_qg),'r-.','Linewidth',1.4); hold on;
loglog(k_qg_low_res,spec_w_qg_low_res,'r:','Linewidth',1.4); hold on;
loglog(k_era5(ind:end),spec_w_era5(ind:end),'k','Linewidth',1.4); hold on;
loglog(k_gcm(indL:end),1.5*k_gcm(indL:end).^(-1/2),'r--','Linewidth',1.4);
loglog(k_gcm(indL:end),0.8*k_gcm(indL:end).^(-4/5),'k--','Linewidth',1.4);

xlabel('Wavenumber k')
ylabel('Amplitude (normalized)')
legend('GCM r=0.01','QG r=0.01','ERA5','k^{-1/2}','k^{-4/5}'); legend boxoff
%legend('GCM','QG r=0.01','QG r=1','QG r=0.01 low res','ERA5','k^{-1/2}','k^{-4/5}','Location','SouthWest'); legend boxoff
legend('GCM','QG r=0.01','QG r=0.01 low res','ERA5','k^{-1/2}','k^{-4/5}','Location','SouthWest'); legend boxoff
%legend('GCM r=0.01','GCM r=0.1','QG r=0.01','QG r=0.01 low res','ERA5','k^{-1/2}','k^{-4/5}'); legend boxoff
title('\rm 1-D w Spectrum');
set(gca,'FontSize',12);
ylim([1e-3 1e1])

% make a plot of the different gcm spectra at different r

load('/disk7/mkohl/inversion/spec_w_gcm_r001.mat');

R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

k_gcm = k*1e6/(2*sqrt(2));

index = [1,4,7,10];

figure

for ii = 1:length(index)
spec_w_gcm = spec_500_r(index(ii),:);
spec_w_gcm = spec_w_gcm/max(spec_w_gcm);

loglog(k_gcm,spec_w_gcm,'Linewidth',1.4); hold on;
end

%loglog(k_qg_r1,spec_w_qg_r1)
xlabel('Wavenumber k')
ylabel('Amplitude')
legend('r=0.01','r=0.1','r=0.4','r=1.0')
legend boxoff
ylim([1e-4 1e1])
set(gca,'FontSize',12)
title('\rm 1-D w spectrum GCM')

clear('k','spec_500_r');