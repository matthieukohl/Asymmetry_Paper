% plot the toy-model results and compare to GCM
% use spec vs lat calculations for the random wavenumber model

close all; clear;

% define r

r = linspace(1e-2,1,300);
%r = linspace(1e-3,1,300);
%r = logspace(-7,0,100);

% define k 

k = [1.8];

% asymmetry

lambda = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
     
     % symmetric rhs    
     [~,lambda(ii,jj)] = toy_model(r(ii),k(jj));
     %lambda(ii,jj) = lambda(ii,jj) + 0.05;
     
     
     % test with a slightly asymmetric rhs 
     %[~,lambda(ii,jj)] = toy_model_asymmetric_rhs(r(ii),k(jj));
     
    end
end

% define gradient colormap
%grad=colorGradient([0.85 0.325 0.098],[0.9290 0.6940 0.1250],size(lambda,2));
%grad=colorGradient([1 0.7812 0.4975],[0 0 0],size(lambda,2));
grad = copper(size(lambda,2));
%grad=colorGradient([0 0 1],[1 0 0],size(lambda,2));

grad = flip(grad,1);


% load GCM fields

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% change interpreter to latex

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

r_gcm = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

% calculate the asymmetry using the full rhs spectrum

% load gcm simulation

lambda_vs_r = zeros(length(r_gcm),1);
skew_rhs = zeros(length(r_gcm),1);

load('/disk7/mkohl/inversion/spec_rhs_gcm_vs_lat.mat');

draws = 100;

for ii = 1:length(r_gcm)
ii
counter = 0;

 for kk = 1:length(lat)
    
spec_rhs_gcm = spec_500_r_vs_lat(ii,kk,:); % select the r=0.01 run
spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);
spec_rhs_gcm = squeeze(spec_rhs_gcm);
spec_rhs_gcm = spec_rhs_gcm';

R = 6371000.0;
Omega = 7.2921e-5;
f = 2 * Omega * sind(lat(kk));
R_eff = cosd(lat(kk))*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];
%Ld = sqrt(S_500_r_vs_lat(ii,kk))*800*1e2/1e-4;
Ld = sqrt(S_500_r_vs_lat(ii,kk))*800*1e2/f;
Ld = Ld/sqrt(2);
Ld = 0.5*Ld; % because our H is the layer height and not the trospopheric height


L = 2*pi*R*cosd(lat(kk));
L = L/Ld;

N = length(spec_rhs_gcm);
x = linspace(0,L,N+1);
dx = L/N;

lambda_pholder = 0;

for jj = 1:draws
    
[rhs] = generate_spectrum(spec_rhs_gcm,L);
skew_rhs(ii) = skew_rhs(ii) + skewness(real(rhs));
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_gcm(ii),real(rhs),dx);


%[~,asymmetry_parameter] = toy_model(r_gcm(ii),3.3);
%
lambda_vs_r(ii) = lambda_vs_r(ii) + asymmetry_parameter;

lambda_pholder = lambda_pholder + asymmetry_parameter;

counter = counter + 1;

end

lambda_vs_lat(kk) = lambda_pholder/draws;

 end

skew_rhs(ii) = skew_rhs(ii)/counter;
lambda_vs_r(ii) = lambda_vs_r(ii)/counter;

 

% sample w-profile

% figure
% plot(x,w_final)
% xlabel('x')
% ylabel('w')
% xlim([x(1) x(end)])
% title(['r=0.01, $\lambda=$',num2str(round(asymmetry_parameter,2))])

end

clear('w_final');


figure(1)

semilogx(r_gcm,lambda_vs_r,'Color','r','Linestyle','--','Linewidth',1.6); hold on;
semilogx(r,lambda,'Color','r','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'Color','b','Linestyle','--','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega,'Color','b','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_mode_omega,'Color','k','Linewidth',1.6); hold on;
legend({'Toy Model Full Spectrum','Toy Model Single Wavenumber','3-D Inversion','GCM','GCM Mode'});
% pos = lgd.Position;
% lgd.FontSize = 10;
% set(lgd,'Position',pos+[0.01 0.02 0 0]);
legend boxoff
set(gca,'FontSize',12)
xlim([r(1) r(end)])
xlabel('Reduction factor r')
ylabel('Asymmetry parameter $\lambda$')
%title('Asymmetry parameter $\lambda$ vs. r')
ylim([0.5 1])
set(gca,'box','off')

%saveas(gcf,'figures_paper/toy_model_asymmetry_vs_gcm','epsc');
