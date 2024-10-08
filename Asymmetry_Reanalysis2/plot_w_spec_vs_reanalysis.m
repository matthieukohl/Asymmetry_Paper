% plot w spec vs reanalysis

close all; clear; 

R = 6371000.0; % Planetary Radius 
dp = 800*1e2;
f = 1e-4;

% ncep2 reanalysis



% load ncep2 data

load('data_theory/theory_ncep2_NH_april29th.mat','r_500_NH','k_NH_500','S_500_NH','spec_NH_500')
load('data_theory/theory_ncep2_SH_april29th.mat','r_500_SH','k_SH_500','S_500_SH','spec_SH_500')

% calculate Ld

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

% calculate k and spec

spec_NH_500_ncep2 = mean(spec_NH_500,2);
spec_SH_500_ncep2 = mean(spec_SH_500,2);

spec_NH_500_ncep2 = spec_NH_500_ncep2/max(spec_NH_500_ncep2);
spec_SH_500_ncep2 = spec_SH_500_ncep2/max(spec_SH_500_ncep2);

N = size(spec_NH_500,1);
L = 2*pi*R*cosd(45);

k = 2*pi/L*[-N/2:N/2-1];
k_NH_ncep2 = k*mean(LD_NH_500);

k = 2*pi/L*[-N/2:N/2-1];
k_SH_ncep2 = k*mean(LD_SH_500);

sum(k_NH_ncep2(N/2:end)'.*spec_NH_500_ncep2(N/2:end))./(sum(spec_NH_500_ncep2(N/2:end)))


clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH','spec_NH_500','spec_SH_500');

% load data erai

load('data_theory/theory_erai_NH_april29th.mat','r_500_NH','k_NH_500','S_500_NH','spec_NH_500')
load('data_theory/theory_erai_SH_april29th.mat','r_500_SH','k_SH_500','S_500_SH','spec_SH_500')

r_500_NH_erai = r_500_NH;
r_500_SH_erai = r_500_SH;

% calculate Ld

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_erai = LD_NH_500;
LD_SH_500_erai = LD_SH_500;

% calculate k and spec

spec_NH_500_erai = mean(spec_NH_500,2);
spec_SH_500_erai = mean(spec_SH_500,2);

spec_NH_500_erai = spec_NH_500_erai/max(spec_NH_500_erai);
spec_SH_500_erai = spec_SH_500_erai/max(spec_SH_500_erai);



N = size(spec_NH_500,1);
L = 2*pi*R*cosd(45);

% LD_NH_500 = 1e6/(2*sqrt(2));
% LD_SH_500 = 1e6/(2*sqrt(2));


k = 2*pi/L*[-N/2:N/2-1];
k_NH_erai = k*mean(LD_NH_500);

k = 2*pi/L*[-N/2:N/2-1];
k_SH_erai = k*mean(LD_SH_500);

k_centroid_erai = sum(k_NH_erai(N/2:end)'.*spec_NH_500_erai(N/2:end))./(sum(spec_NH_500_erai(N/2:end)))



clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH','spec_NH_500','spec_SH_500');

% load data era5

% native resolution

% load('data_theory/theory_NH_dec7th.mat','r_500_NH','k_NH_500','S_500_NH','spec_NH_500')
% load('data_theory/theory_SH_dec7th.mat','r_500_SH','k_SH_500','S_500_SH','spec_SH_500')

% coarse grained

load('data_theory/theory_era5_NH_coarse_grained_july_5th.mat','r_500_NH','k_NH_500','S_500_NH','spec_NH_500')
load('data_theory/theory_era5_SH_coarse_grained_july_5th.mat','r_500_SH','k_SH_500','S_500_SH','spec_SH_500')


r_500_NH_era5 = r_500_NH;
r_500_SH_era5 = r_500_SH;



LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

LD_NH_500_era5 = LD_NH_500;
LD_SH_500_era5 = LD_SH_500;

% calculate k and spec

spec_NH_500_era5 = mean(spec_NH_500,2);
spec_SH_500_era5 = mean(spec_SH_500,2);

spec_NH_500_era5 = spec_NH_500_era5/max(spec_NH_500_era5);
spec_SH_500_era5 = spec_SH_500_era5/max(spec_SH_500_era5);

N = size(spec_NH_500,1);
L = 2*pi*R*cosd(45);

% LD_NH_500 = 1e6/(2*sqrt(2));
% LD_SH_500 = 1e6/(2*sqrt(2));

k = 2*pi/L*[-N/2:N/2-1];
k_NH_era5 = k*mean(LD_NH_500);

k = 2*pi/L*[-N/2:N/2-1];
k_SH_era5 = k*mean(LD_SH_500);

k_centroid_era5 = sum(k_NH_era5(N/2:end)'.*spec_NH_500_era5(N/2:end))./(sum(spec_NH_500_era5(N/2:end)))

clear('r_500_NH','r_500_SH','k_500_NH','k_500_SH','S_500_NH','S_500_SH','spec_NH_500','spec_SH_500');


% make a plot

semilogx(k_NH_ncep2,spec_NH_500_ncep2,'g'); hold on;
semilogx(k_NH_erai,spec_NH_500_erai,'r'); hold on;
semilogx(k_NH_era5,spec_NH_500_era5,'b');
xlabel('Reduction factor r')
ylabel('Spectrum')
legend('ncep2','erai','era5'); legend boxoff

figure
loglog(k_NH_ncep2,spec_NH_500_ncep2,'g'); hold on;
loglog(k_NH_erai,spec_NH_500_erai,'r'); hold on;
loglog(k_NH_era5,spec_NH_500_era5,'b'); hold on;
loglog(k_NH_era5,1e2*k_NH_era5.^(-3),'k--');
xlabel('Reduction factor r')
ylabel('Spectrum')
ylim([1e-10 1e1])
legend('ncep2','erai','era5'); legend boxoff


