% plot the toy-model results and compared to GCM

close all; clear;

% define r

r = linspace(1e-2,1,300);
%r = linspace(1e-3,1,300);
%r = logspace(-7,0,100);

% define k 

%k = 1:0.3:1.9;
%k = 1.3:0.3:2.2;

%k = 1.2:0.3:2.1;
%k = [0.73,1.5,1.8];
%k = [0.73,1.5,2.3];
k = [1.2,1.5,1.7];
k = [k,6.6];

%k = [1.8,2.3,6.6,11];


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

load('/disk7/mkohl/inversion/spec_rhs_gcm_r001.mat');

Ld = 1e6/(2*sqrt(2));
%Ld = 5.79*1e5; % using S500 from model
k_gcm = k*Ld;

draws = 100;

for ii = 1:length(r_gcm)
    
spec_rhs_gcm = spec_500_r(ii,:); % select the r=0.01 run
spec_rhs_gcm = spec_rhs_gcm/max(spec_rhs_gcm);

R = 6371000.0;
%lat = 41.2336;
lat = 45;
L = 2*pi*R*cosd(lat);
L = L/Ld;

N = length(spec_rhs_gcm);
x = linspace(0,L,N+1);
dx = L/N;

for jj = 1:draws
    
[rhs] = generate_spectrum(spec_rhs_gcm,L);
skew_rhs(ii) = skew_rhs(ii) + skewness(real(rhs));
[w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(r_gcm(ii),real(rhs),dx);


%[~,asymmetry_parameter] = toy_model(r_gcm(ii),3.3);
%
lambda_vs_r(ii) = lambda_vs_r(ii) + asymmetry_parameter;



end

skew_rhs(ii) = skew_rhs(ii)/draws;
lambda_vs_r(ii) = lambda_vs_r(ii)/draws;



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
x0=10; y0=10; width=1300; height=300; set(gcf,'position',[x0,y0,width,height])

%x0=10; y0=10; width=400; height=300; set(gcf,'position',[x0,y0,width,height])

grayColor = [.7 .7 .7];

subplot(1,3,2)
h = semilogx(r,lambda,'Linewidth',1.6); hold on;
%semilogx(r_gcm,lambda_turb_omega,'Color',grayColor,'Linewidth',1.2); hold on;
%semilogx(r_gcm,lambda_turb_omega_QG,'Color',grayColor,'Linestyle','--','Linewidth',1.2); hold on;
%semilogx(r_gcm,lambda_vs_r,'Color',grayColor,'Linestyle',':','Linewidth',1.2);
set(h,{'color'},num2cell(grad,2)); % apply gradient colormap
%legend('k=1.0','k=1.3','k=1.6','k=1.9','k=6.0');
%legend('k=1.2','k=1.5','k=1.8','k=2.1','k=6.0');
%lgd = legend({'k=0.73','k=1.5','k=2.3','k=6.6','GCM','3-D Inversion'});
lgd = legend({'k=1.2','k=1.5','k=1.7','k=6.6'},'Location','NorthEast');
%legend('k=1.8','k=2.3','k=6.6','k=11');
%legend('k=1.2','k=1.5','k=1.8');
%pos = lgd.Position;
%lgd.FontSize = 10;
%set(lgd,'Position',pos+[0.01 0.02 0 0]);
legend boxoff
set(gca,'FontSize',12)
%xticks([0.001 0.01 0.1 1])
%xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^{0}'})
xlim([r(1) r(end)])
xlabel('Reduction factor r')
ylabel('Asymmetry parameter $\lambda$')
title('Asymmetry parameter $\lambda$ vs. r')
ylim([0.5 1])
set(gca,'box','off')

%saveas(gcf,'figures_paper/toy_model_asymmetry_and_w_profile_alternative_one_panel','epsc');

% calculate the vertical velocity profile at fixed r

k = [1.2,1.5,1.7];
k = [k,6.6];

%[~,k1p6] = min(abs(k-1.6));
[~,k1p6] = min(abs(k-1.7));
%[~,k1p6] = min(abs(k-0.73));

r_sample = [0.01,0.1,0.5,1.0];

w_profile = zeros(300,length(r_sample));

for ii = 1:length(r_sample)
[w_profile(:,ii),~] = toy_model(r_sample(ii),k(k1p6));
end

x = linspace(0,2*pi/k(k1p6),300);

% define gradient colormap
%grad=colorGradient([0.85 0.325 0.098],[0.9290 0.6940 0.1250],size(w_profile,2));
grad = copper(size(w_profile,2));
grad = flip(grad,1);

subplot(1,3,1)
h = plot(x,w_profile,'Linewidth',1.6); 
set(h,{'color'},num2cell(grad,2)); % apply gradient colormap
legend('r=0.01','r=0.1','r=0.5','r=1.0','Location','NorthWest');
legend boxoff
set(gca,'FontSize',12)
xlabel('x')
ylabel('Vertical velocity w')
title('w profile (k=1.7)')
%title('w profile (k=0.73)')
set(gca,'box','off')
xlim([x(1) x(end)])

% plot lambda vs k for fixed r

% asymmetry
r = [0.01,0.1,0.5];
k = linspace(1,11,20);

lambda = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
     [~,lambda(ii,jj)] = toy_model(r(ii),k(jj));
        
    end
end

[a,ind] = min(abs(r-0.01));

subplot(1,3,3)

grad = copper(length(r));
grad = flip(grad,1);

h = plot(k,lambda,'Linewidth',1.4);
set(h,{'color'},num2cell(grad,2));
xlabel('Wavenumber k')
ylabel('Asymmetry parameter $\lambda$');
lgd = legend('r=0.01','r=0.1','r=0.5','Location','NorthWest'); 
legend boxoff;
lgd.FontSize = 11;
xlim([1 11]);
ylim([0.5 1])
set(gca,'box','off')
set(gca,'FontSize',12);
title('Asymmetry parameter $\lambda$ vs. k')

annotation('textbox', [0.09, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.37, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')
annotation('textbox', [0.65, 0.94, 0.01, 0.03], 'String', "(c)",'EdgeColor','none')

set(groot, 'defaultTextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%saveas(gcf,'/disk7/mkohl/figures_paper/toy_model_asymmetry_and_w_profile_three_panels_with_inversion','epsc');
%saveas(gcf,'/disk7/mkohl/figures_paper/toy_model_asymmetry_and_w_profile_three_panels_with_inversion.pdf');
%saveas(gcf,'/net/graupel/archive1/users/mkohl/inversion/figures_paper/toy_model_asymmetry_and_w_profile_three_panels_with_inversion','epsc');
saveas(gcf,'figures_paper/toy_model_asymmetry_and_w_profile_three_panels','epsc');
%saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_abs_rhs','epsc');


% %test
% 
% lambda = zeros(length(r),100);
% 
% for ii = 1:length(r)
%     for jj = 1:100
%         
%      % test with a slightly asymmetric rhs 
%      [~,lambda(ii,jj)] = toy_model_asymmetric_rhs(r(ii),2);
%      
%     end
% end
