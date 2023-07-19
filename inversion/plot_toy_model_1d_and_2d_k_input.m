% plot the 1d and 2d sine and box model where
% we use ku from the centroid of the w spectrum as input

close all; clear;

%r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];
r = logspace(-2,0,10);

k = [1.8];
k2 = k^2;

% box theory

lambda_1d = zeros(length(r),1);
lambda_2d = zeros(length(r),1);

for ii = 1:length(r)
        
% 1d 

lambda_1d(ii) = toy_model_1d_non_sine_refined(r(ii),k);

% 2d

lambda_2d(ii) = toy_model_2d_non_sine_refined(r(ii),k);
    
end



% sine theory


lambda_1d_sin = zeros(length(r),1);
lambda_2d_sin = zeros(length(r),1);



for ii = 1:length(r)
 
        
% 1d 

[a,lambda_1d_sin(ii)] = toy_model(r(ii),k);

% 2d

[a,lambda_2d_sin(ii)] = toy_model_2d_alternative(r(ii),k);

end


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


set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

figure
semilogx(r,lambda_2d_sin,'r','Linewidth',1.4); hold on;
semilogx(r,lambda_2d,'r--','Linewidth',1.4); hold on;
semilogx(r,lambda_1d_sin,'b','Linewidth',1.4); hold on;
semilogx(r,lambda_1d,'b--','Linewidth',1.4); hold on;

%semilogx(r_gcm,lambda_turb_omega,'k'); hold on;
%semilogx(r_gcm,lambda_turb_omega_QG,'k--')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter $\lambda$')
legend('2d sine','2d box','1d sine','1d box','Location','SouthWest'); legend boxoff;
set(gca,'FontSize',12)
set(gca,'box','off')
ylim([0.5 0.9])

saveas(gcf,'figures_paper/1d_vs_2d_box_and_sine_model','epsc');