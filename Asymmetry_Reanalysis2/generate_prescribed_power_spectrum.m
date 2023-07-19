% generate a prescribed power spectrum

close all; clear;

% load gcm lambda values

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% change interpreter to latex

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex')

% load r simulation results

load(['/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_r_',low,'_',averaging,'.mat'])


% load gcm simulation

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];
lambda_vs_r = zeros(length(r),1);
lambda = zeros(length(r),1);
lambda_toy = zeros(length(r),1);
skew_rhs = zeros(length(r),1);

load('/disk7/mkohl/inversion/spec_rhs_gcm_r001.mat');

Ld = 1e6/(2*sqrt(2));
%Ld = 5.79*1e5; % using S500 from model
k_gcm = k*Ld;

draws = 20;

for ii = 1:length(r)
    
spec_rhs_gcm = spec_500_r(ii,:); % select the r=0.01 run
spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);

R = 6371000.0;
lat = 41.2336;
L = 2*pi*R*cosd(lat);
L = L/Ld;

N = length(spec_rhs_gcm);
x = linspace(0,L,N+1);
dx = L/N;

for jj = 1:draws
    
[rhs] = generate_spectrum(spec_rhs_gcm,L);
skew_rhs(ii) = skew_rhs(ii) + skewness(real(rhs));
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r(ii),real(rhs),dx);

%
lambda_vs_r(ii) = lambda_vs_r(ii) + asymmetry_parameter;
lambda(ii) = lambda(ii) + asymmetry(w_final(5:end-5));
end

skew_rhs(ii) = skew_rhs(ii)/draws;
lambda_vs_r(ii) = lambda_vs_r(ii)/draws;
lambda(ii) = lambda(ii)/draws;

lambda_toy(ii) = toy_model_nondimensional(r(ii),1.8);

% sample w-profile

% figure
% plot(x,w_final)
% xlabel('x')
% ylabel('w')
% xlim([x(1) x(end)])
% title(['r=0.01, $\lambda=$',num2str(round(asymmetry_parameter,2))])

end

figure
semilogx(r,skew_rhs)
xlabel('Skew rhs')
ylabel('\lambda')
set(gca,'FontSize',12);

figure
semilogx(r,lambda_turb_omega,'b'); hold on;
semilogx(r,lambda_turb_omega_QG,'b--');
semilogx(r,lambda,'r'); hold on;
semilogx(r,lambda_toy,'r--');
legend('GCM','3D Inversion','Theory Full Spectrum','Theory Single k');
legend boxoff
xlabel('Reduction factor r')
ylabel('$\lambda$')
set(gca,'FontSize',12);
title('r-simulations')


% test it on QG simulations

load('spec_rhs_qg_r001.mat')

k_qg = kx;
spec_rhs_qg = power_rhs_kx;

spec_rhs_qg = spec_rhs_qg/max(spec_rhs_qg);

[a,ind_qg] = min(abs(kx-36));

r = 0.01;
N = length(spec_rhs_qg);
L = 12*pi;
dx = L/N;
x = linspace(0,L,N+1);

clear('kx','power_rhs_kx');

lambda_qg = 0;

for ii = 1:draws

[rhs] = generate_spectrum(spec_rhs_qg,L);
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r,real(rhs),dx);

lambda_qg = lambda_qg + asymmetry_parameter;
end

lambda_qg = lambda_qg/draws;

figure
plot(x,w_final)
xlabel('x')
ylabel('w')
xlim([x(1) x(end)])
title('QG r=0.01, $\lambda=0.91$')

% now do the ERA5 comparison

% era5

load('data_era5/asymmetry_era5_NH_x.mat'); 
load('data_era5/asymmetry_era5_NH_xt.mat'); 

lv = 500;

latN = 70; latS = 30;

[~,lat30] = min(abs(lat-latS)); 
[~,lat70] = min(abs(lat-latN)); 
 
[~,p500] = min(abs(level-lv));

lambda_era5_NH_x = squeeze(mean(lambda_era5_NH_x(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_era5/asymmetry_era5_SH_x.mat'); 
load('data_era5/asymmetry_era5_SH_xt.mat');

[~,latm30] = min(abs(lat+latS));
[~,latm70] = min(abs(lat+latN));

[~,p500] = min(abs(level-lv));

lambda_era5_SH_x = squeeze(mean(lambda_era5_SH_x(latm30:latm70,p500,:),1));
lambda_era5_SH_xt = squeeze(mean(lambda_era5_SH_xt(latm30:latm70,p500,:),1));

clear('level','lat','p500','lat30','lat70')

load('data_theory/theory_NH_dec7th.mat');
load('data_theory/theory_SH_dec7th.mat');

dp = 800*1e2;
f = 1e-4;

LD_NH_500 = sqrt(S_500_NH)*dp/f;
LD_NH_500 = LD_NH_500/sqrt(2); % due to our choice of nondimensionalization
LD_NH_500 = 0.5*LD_NH_500; % because our H is the layer and not the tropospheric height

LD_SH_500 = sqrt(S_500_SH)*dp/f; 
LD_SH_500 = LD_SH_500/sqrt(2); % due to our choice of nondimensionalization
LD_SH_500 = 0.5*LD_SH_500; % because our H is the layer not the tropospheric height

load('spec_rhs_era5_new.mat');
spec_rhs_era5_NH = spec_rhs_era5;
k_era5 = k;

clear('spec_rhs_era5','k');


load('spec_rhs_era5_by_month_sh.mat');
spec_rhs_era5_SH = spec_rhs_era5;


% spec_rhs_era5 = spec_rhs_winter;
% 
% spec_rhs_era5 = spec_rhs_era5/max(spec_rhs_era5);


lambda_era5_NH = zeros(12,1);
lambda_era5_SH = zeros(12,1);

for ii = 1:12
    
R = 6371000.0;
lat = 41.2336;
L = 2*pi*R*cosd(lat);
%Ld = 1e6/(2*sqrt(2)); % nondimensionalize
Ld = LD_NH_500(ii); % use Ld with S500 from model
L = L/Ld;

%N = length(spec_rhs_era5);
N = size(spec_rhs_era5_NH,2);
x = linspace(0,L,N+1);
dx = L/N;
    

 for jj = 1:draws

[rhs] = generate_spectrum(mean(spec_rhs_era5_NH,1),L);
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_500_NH(ii),real(rhs),dx);

lambda_era5_NH(ii) = lambda_era5_NH(ii) + asymmetry_parameter;

[rhs] = generate_spectrum(mean(spec_rhs_era5_SH,1),L);

[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_500_SH(ii),real(rhs),dx);

lambda_era5_SH(ii) = lambda_era5_SH(ii) + asymmetry_parameter;


 end
 
lambda_era5_NH(ii) = lambda_era5_NH(ii)/draws;
lambda_era5_SH(ii) = lambda_era5_SH(ii)/draws;

end

lambda_era5_NH_theory = lambda_era5_NH;
lambda_era5_SH_theory = lambda_era5_SH;


hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b'); hold on;
plot(lambda_era5_NH_theory,'k--','linewidth',1.5); hold on;
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2','ERAI','ERA5','Location','Northwest','Turb','Modal'); legend boxoff
l1 = legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[-0.005 0 0 0]);
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','PG','Location','Northwest'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

plot(lambda_era5_SH_x,'linewidth',1.5,'Color','b'); hold on;
%plot(lambda_era5_SH_xt,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
plot(lambda_era5_SH_theory,'k--','linewidth',1.5); hold on;
%plot(lambda_SH_modal_ZG,'linewidth',1.5,'Color','k','Linestyle','--');hold on
%plot(lambda_SH_modal_PG,'linewidth',1.5,'Color','k','Linestyle','-.');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
legend('ERA5','Location','Northwest','1-D Toy Model','1-D Modal Theory'); legend boxoff
title(['$ \lambda$ (',num2str(lv),'hPa): ',num2str(latS),'$^{\circ}$-',num2str(latN),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

%figure
% plot(x,w_final)
% xlabel('x')
% ylabel('w')
% xlim([x(1) x(end)])
% title('ERA5, $\lambda=0.63$')
% 
% 
% N = length(k_gcm);
% dx = L/N;
% P = spec_rhs_gcm;
% P = P';
% 
% % L = 4;
% % N = 1000;
% % dx = L/N;
% % k_gcm = 2*pi/L*[0:N/2,-N/2+1:-1];
% % P = zeros(N,1);
% % P(3) = 1; P(end-1) = 1;
% 
% 
% A = sqrt(P/(2*N*dx)) .* exp( 1i * rand(N,1) );
% 
% 
% 
% %N = 1e3;
% %P = ones(N,1);
% %k1 = 2*pi*[0:N/2,-N/2+1:-1];
% %P = zeros(N,1);
% %P(3) = 1; P(end-1) = 1;
% %P(k1>0) = k1(k1>0); P(k1<0) = -k1(k1<0); %P=P';
% %Ts = 1/1000;
% %A = sqrt(P/(2*N*Ts)) .* exp( 1i * rand(N,1) );
% %A = sqrt(P/(2*N*Ts)) .* exp( 2*pi*1i * rand(N,1) );
% 
% 
% % DC sample must be real.
% n0 = floor( N/2 ) + 1;
% A(n0) = P(n0);
% 
% % For even-length cases, the first sample has no
% % equivalent positive frequency, so we force it
% % to be real too.
% A(1) = abs( A(1) );
% 
% % Force the spectrum to be conjugate symetric.
% if( mod(N,2) == 0 )
%   A(2:n0-1) = flipud( conj(A(n0+1:end)) );
% else
%   A(1:n0-1) = flipud( conj(A(n0+1:end)) );
% end
% 
% % Time domain data.
% x = ifft(ifftshift(A));
% 
% % % test
% % k = 2*pi*(-N/2:N/2-1);
% % power = 2*abs(fft(x)).^2;
% % 
% % %plot(k,k,'k'); hold on;
% % %loglog(k,power,'b'); hold on;
% % %loglog(k1,P,'k--');
% 
% % test
% k = 2*pi/L*(-N/2:N/2-1);
% power = 2*L*abs(fft(x)).^2;
% %power = power/max(power);
% 
% %plot(k,k,'k'); hold on;
% loglog(k,power,'b'); hold on;
% loglog(k_gcm,P,'k--');
% 
% 
% 
% 
% % test arbsnoise
% 
% L = 2*pi;
% N = 1000;
% k = 2*pi/L*[0:N/2,-N/2+1:-1];
% alpha = 2;
% 
% counter = 0;
% for ii = 1:50
% f = arbssnoise(N,alpha);
% power = power + abs(fft(f)).^2;
% counter = counter + 1;
% end
% power = power/counter;
% 
% loglog(k,100*power); hold on;
% loglog(k,k.^2);
% 
% % code
% N = 1e3;
% P = ones(N,1);
% Ts = 1/1000;
% A = sqrt(P/(2*N*Ts)) .* exp( 1j * rand(N,1) );
% 
% % % DC sample must be real.
% n0 = floor( N/2 ) + 1;
% %A(n0) = P(n0);
% A(n0) = abs(A(n0));
% % 
% % For even-length cases, the first sample has no
% % equivalent positive frequency, so we force it
% % to be real too.
% A(1) = abs( A(1) );
% 
% % Force the spectrum to be conjugate symetric.
% if( mod(N,2) == 0 )
%   A(2:n0-1) = flipud( conj(A(n0+1:end)) );
% else
%   A(1:n0-1) = flipud( conj(A(n0+1:end)) );
% end
% 
% % Time domain data.
% x = ifft( ifftshift( A ) );
% 
% 
% k = 2*pi*[0:N/2,-N/2+1:-1];
% 
% plot(k,abs(fft(x)).^2);
