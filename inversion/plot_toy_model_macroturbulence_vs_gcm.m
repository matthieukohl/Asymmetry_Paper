% plot the toy-model results and compared to GCM

close all; clear;

% define r

%r = linspace(1e-2,1,300);
%r = linspace(1e-3,1,300);
%r = logspace(-7,0,100);

% define k 

%k = [1.8];
%k = [2.8];

% calculate k

load('spec_rhs_gcm_r001.mat');
load('figures_paper/deformation_radius_r_simulation.mat');
[a,p500] = min(abs(level-500e2));

r = r_gcm;
 
% Ld = squeeze(Ld_vs_r(p500,:));
% N = length(k);
% 
% 
% for tt = 1:10
% k_gcm = k(1:N/2);   
% k_gcm = k_gcm*Ld(tt);
% %k_gcm = k_gcm*1e6/(2*sqrt(2));
% 
% k_test(tt) = sum(k_gcm.*spec_500_r(tt,1:N/2))./(sum(spec_500_r(tt,1:N/2)));
% 
% end



% % asymmetry
% 
% lambda = zeros(length(r),length(k));
% 
% for ii = 1:length(r)
%     for jj = 1:length(k)
%      
%      % symmetric rhs    
%      [~,lambda(ii,jj)] = toy_model(r(ii),k(jj));
%      %lambda(ii,jj) = lambda(ii,jj) + 0.05;
%      
%      
%      % test with a slightly asymmetric rhs 
%      %[~,lambda(ii,jj)] = toy_model_asymmetric_rhs(r(ii),k(jj));
%      
%     end
% end

% asymmetry

lambda = zeros(length(r),1);
lambda_k_rhs = zeros(length(r),1);

load('figures_paper/k_w_spectrum_r_simulation_gcm.mat');
k_w = k_500_r.*Ld;
k_w = round(k_w,1);
clear('k_500_r','Ld');

load('figures_paper/k_rhs_spectrum_r_simulation_gcm.mat');
k_rhs = k_500_r.*Ld;
k_rhs = round(k_rhs,1);



for ii = 1:length(r)
     
     % symmetric rhs    
     %[~,lambda(ii)] = toy_model(r(ii),k_test(ii));
     %[~,lambda(ii)] = toy_model(r(ii),1.7);
     %[~,lambda_k_rhs(ii)] = toy_model(r(ii),3.2);
     
     [~,lambda(ii)] = toy_model(r(ii),k_w(ii));
     [~,lambda_k_rhs(ii)] = toy_model(r(ii),k_rhs(ii));
     
     %[~,lambda(ii)] = toy_model_2d_alternative(r(ii),k_test(ii));
     
     %lambda(ii,jj) = lambda(ii,jj) + 0.05;
     
     
     % test with a slightly asymmetric rhs 
     %[~,lambda(ii,jj)] = toy_model_asymmetric_rhs(r(ii),k(jj));
     
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

%load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])
load(['figures_paper/lambda_r_',low,'_',averaging,'_rms1em3_relative_error.mat'])


r_gcm = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

% % calculate the asymmetry using the full rhs spectrum
% 
% % load gcm simulation
% 
% lambda_vs_r = zeros(length(r_gcm),1);
% skew_rhs = zeros(length(r_gcm),1);
% 
% load('/disk7/mkohl/inversion/spec_rhs_gcm_r001.mat');
% load('figures_paper/deformation_radius_r_simulation_p700_p400_average.mat');
% 
% %Ld = 1e6/(2*sqrt(2));
% %Ld = 5.79*1e5; % using S500 from model; see file plot_w_spec_r_simulation
% %Ld = Ld*2*sqrt(2)/sqrt(8); % define using parabolic w-profile
% %Ld = 5.21*1e5; % using weighted-sine method
% %k_gcm = k*Ld;
% 
% draws = 1000;
% 
% for ii = 1:length(r_gcm)
%     
% ii
%     
% [pholder,p500] = min(abs(level-50000));
% Ld = Ld_vs_r(p500,ii);
% %Ld = 1e6/(2*sqrt(2));
% %Ld = 3.5*1e5;
% 
% 
% k_gcm = k*Ld;
%     
% spec_rhs_gcm = spec_500_r(ii,:); % select the r=0.01 run
% spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);
% 
% R = 6371000.0;
% %lat = 41.2336;
% lat = 45.5; % gcm domain doesn't start/end exactly at 25-65 lat;
% L = 2*pi*R*cosd(lat);
% L = L/Ld;
% 
% N = length(spec_rhs_gcm);
% x = linspace(0,L,N+1);
% dx = L/N;
% 
% w_aggregate = [];
% for jj = 1:draws
%         
% [rhs] = generate_spectrum(spec_rhs_gcm,L);
% skew_rhs(ii) = skew_rhs(ii) + skewness(real(rhs));
% [w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_gcm(ii),real(rhs),dx);
% 
% 
% %[~,asymmetry_parameter] = toy_model(r_gcm(ii),3.3);
% %
% lambda_vs_r(ii) = lambda_vs_r(ii) + asymmetry_parameter;
% 
% w_aggregate = cat(1,w_aggregate,w_final);
% 
% end
% 
% skew_rhs(ii) = skew_rhs(ii)/draws;
% lambda_vs_r(ii) = lambda_vs_r(ii)/draws;
% 
% 
% 
% % sample w-profile
% 
% % figure
% % plot(x,w_final)
% % xlabel('x')
% % ylabel('w')
% % xlim([x(1) x(end)])
% % title(['r=0.01, $\lambda=$',num2str(round(asymmetry_parameter,2))])
% 
% end
% 
% clear('w_final');

set(groot,'DefaultTextInterpreter','default');
set(groot,'DefaultTextInterpreter','default');


figure(1)

%semilogx(r_gcm,lambda_vs_r,'Color','r','Linestyle','--','Linewidth',1.6); hold on;
semilogx(r,lambda,'Color','r','Linewidth',1.6); hold on;
semilogx(r,lambda_k_rhs,'Color','r','Linestyle',':','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'Color','b','Linestyle','--','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega,'Color','b','Linewidth',1.6); hold on;
%semilogx(r_gcm,lambda_mode_omega,'Color','k','Linewidth',1.6); hold on;
%legend({'Toy Model Full Spectrum','Toy Model Single Wavenumber','3-D Inversion','GCM','GCM Mode'});
%legend({'1-D Toy Model ($k_{w}$)','1-D Toy Model ($k_{RHS}$)','3-D Inversion','GCM'},'FontSize',12,'Interpreter','latex');
%legend({'1-D Toy Model (k=1.7)','1-D Toy Model (k=3.2)','3-D Inversion','GCM'},'FontSize',12,'Interpreter','latex');
legend({'1-D Toy Model (k from w-spectrum)','1-D Toy Model (k from Adv-spectrum)','3-D Inversion','GCM'},'FontSize',12,'Interpreter','latex');
% pos = lgd.Position;
% lgd.FontSize = 10;
% set(lgd,'Position',pos+[0.01 0.02 0 0]);
legend boxoff
set(gca,'FontSize',12)
xlim([r(1) r(end)])
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
%title('Asymmetry parameter $\lambda$ vs. r')
ylim([0.5 1])
set(gca,'box','off')

%saveas(gcf,'figures_paper/toy_model_asymmetry_vs_gcm_kw_and_krhs','epsc');
saveas(gcf,'/net/graupel/archive1/users/mkohl/figures_paper_asymmetry/toy_model_asymmetry_vs_gcm_kw_and_krhs_continous_r','epsc');
%saveas(gcf,'/net/graupel/archive1/users/mkohl/figures_paper_asymmetry/toy_model_asymmetry_vs_gcm_kw_and_krhs','pdf');
