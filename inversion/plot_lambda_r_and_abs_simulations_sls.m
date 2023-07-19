% plot lambda for w and w_qg for r and abs simulations

clear; close all;

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% change interpreter to latex

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');


% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

lambda_mode_omega_QG_rhs(end) = 0.5;

figure 

h1 = semilogx(r,lambda_turb_omega,'b','linewidth',1.2); hold on;
h2 = semilogx(r,lambda_turb_omega_QG,'b--','linewidth',1.2); hold on;
h3 = semilogx(r,lambda_mode_omega,'r','linewidth',1.2); hold on;
h4 = semilogx(r,lambda_mode_omega_QG,'r--','linewidth',1.2); hold on;
x1 = xlabel('Reduction factor r');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
y1 = ylabel('Asymmetry parameter $\lambda$');
%y1 = ylabel(sprintf('Asymmetry parameter \\lambda'));
%y1pos = get(y1,'position');
%l1 = legend({'GCM_{turbul.}','QG_{turbul.}','GCM_{mode}','QG_{mode}'},'Location','NorthEast','FontSize',8.5); legend boxoff
%l1 = legend([h3,h4,h1,h2],{'GCM$_{turbul.}$','QG$_{turbul.}$','GCM$_{mode}$','QG$_{mode}$'},'Location','NorthEast','FontSize',10); legend boxoff
%l1 = legend([h3,h4,h1,h2],{'GCM$_{mode}$','QG$_{mode}$','GCM$_{turbul.}$','QG$_{turbul.}$'},'Location','NorthEast','FontSize',10); legend boxoff
%l1pos = get(l1,'position');
%set(l1,'position',l1pos+[0.02 0.01 0 0]);
%t1 = title('Inversion with r($\omega$)');
%t1 = title('\rm Inversion with r(w)');
%t1 = title('$\nabla^2(r(w)N^2(z)w)+f_0^2w_{zz} = RHS$','interpreter','latex');
%t1pos = get(t1,'position');
%set(t1,'position',t1pos+[0 0.01 0]);
set(gca,'box','off')
ylim([0.5 1])
%set(y1,'position',y1pos+[-0.0023 0.05 0]);
set(gca,'FontSize',11);
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_r','epsc');

figure 

semilogx(r,lambda_turb_omega,'b','linewidth',1.2); hold on;
semilogx(r,lambda_turb_omega_QG_rhs,'b--','linewidth',1.2); hold on
semilogx(r,lambda_mode_omega,'r','linewidth',1.2); hold on;
semilogx(r,lambda_mode_omega_QG_rhs,'r--','linewidth',1.2);
x1 = xlabel('Reduction factor r');
y2 = ylabel('Asymmetry parameter $\lambda$');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
%ylabel('\lambda')
set(gca,'FontSize',11);
%t2 = title('Inversion with r($\omega$)=1');
%t2 = title('\rm Inversion with r(w)=1');
%set(t2,'position',t2pos+[0 0.01 0]);
set(gca,'box','off')
ylim([0.5 1])
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_r_rhs','epsc');


% load abs results

load(['figures_paper/lambda_abs_',low,'_',averaging,'.mat'])

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

figure

plot(surf_temp,lambda_turb_omega,'b','linewidth',1.2); hold on;
plot(surf_temp,lambda_turb_omega_QG,'b--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega,'r','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega_QG,'r--');
x1 = xlabel('Global-mean surface air temperature (K)');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
y2 = ylabel('Asymmetry parameter $\lambda$');
%y2 = ylabel(sprintf('Asymmetry parameter \\lambda %g\t'));
%y2pos = get(y2,'position');
%set(y2,'position',y2pos+[-5 0 0]);
set(gca,'FontSize',11);
%title('\rm RHS+r') 
ylim([0.5 1])
set(gca,'box','off')
ylim([0.5 1])
L(1) = plot(nan, nan, 'k');
L(2) = plot(nan, nan, 'k--');
legend(L, {'Simulation', 'Moist Omega Equation'},'box','off','FontSize',14,'Location','NorthWest')
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_abs','epsc');



figure

plot(surf_temp,lambda_turb_omega,'b','linewidth',1.2); hold on;
plot(surf_temp,lambda_turb_omega_QG_rhs,'b--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega,'r','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega_QG_rhs,'r--','linewidth',1.2);
x2 = xlabel('Global-mean surface air temperature (K)');
y2 = ylabel('Asymmetry parameter $\lambda$');
%x2pos = get(x2,'position');
%set(x2,'position',x2pos+[0 -0.03 0]);
%set(gca,'box','off')
set(gca,'box','off')
%ylabel('\lambda');
set(gca,'FontSize',11);
%title('\rm RHS') 
ylim([0.5 1])
L(1) = plot(nan, nan, 'k');
L(2) = plot(nan, nan, 'k--');
legend(L, {'Simulation', 'Moist Omega Equation with r=1'},'box','off','FontSize',14,'Location','NorthWest')
lgd.TextColor = 'black';
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_abs_rhs','epsc');


set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex')