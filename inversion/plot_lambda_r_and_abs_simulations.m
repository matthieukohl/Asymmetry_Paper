% plot lambda for w and w_qg for r and abs simulations

clear; close all;

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% change interpreter to latex

set(groot, 'DefaultAxesTickLabelInterpreter','latex');
set(groot, 'DefaultLegendInterpreter','latex');


% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

%load(['figures_paper/lambda_r_',low,'_',averaging,'_rms1em5.mat'])


r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

figure(1)
x0=10; y0=10; width=950; height=600; set(gcf,'position',[x0,y0,width,height])

subplot(2,2,1)
h1 = semilogx(r,lambda_turb_omega,'b','linewidth',1.2); hold on;
h2 = semilogx(r,lambda_turb_omega_QG,'b--','linewidth',1.2); hold on;
h3 = semilogx(r,lambda_mode_omega,'r','linewidth',1.2); hold on;
h4 = semilogx(r,lambda_mode_omega_QG,'r--','linewidth',1.2); hold on;
x1 = xlabel('Reduction factor r');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
y1 = ylabel('Asymmetry parameter $\lambda$');
%y1 = ylabel(sprintf('Asymmetry parameter \\lambda'));
y1pos = get(y1,'position');
%l1 = legend({'GCM_{turbul.}','QG_{turbul.}','GCM_{mode}','QG_{mode}'},'Location','NorthEast','FontSize',8.5); legend boxoff
%l1 = legend([h3,h4,h1,h2],{'GCM$_{mode}$','QG$_{mode}$','GCM$_{turbul.}$','QG$_{turbul.}$'},'Location','NorthEast','FontSize',10); legend boxoff
l1 = legend([h3,h4,h1,h2],{'Mode','Inverted Mode','Turb.','Inverted Turb.'},'Location','NorthEast','FontSize',10); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[0.02 0.01 0 0]);
t1 = title('QG Inversion uses r($\omega$)');
%t1 = title('\rm Inversion with r(w)');
%t1 = title('$\nabla^2(r(w)N^2(z)w)+f_0^2w_{zz} = RHS$','interpreter','latex');
t1pos = get(t1,'position');
set(t1,'position',t1pos+[0 0.01 0]);
set(gca,'box','off')
ylim([0.5 1])
set(y1,'position',y1pos+[-0.0023 0.05 0]);
set(gca,'FontSize',11);

subplot(2,2,2)
semilogx(r,lambda_turb_omega,'b','linewidth',1.2); hold on;
semilogx(r,lambda_turb_omega_QG_rhs,'b--','linewidth',1.2); hold on
semilogx(r,lambda_mode_omega,'r','linewidth',1.2); hold on;
semilogx(r,lambda_mode_omega_QG_rhs,'r--','linewidth',1.2);
x1 = xlabel('Reduction factor r');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
%ylabel('\lambda')
set(gca,'FontSize',11);
t2 = title('QG Inversion uses r($\omega$)=1');
%t2 = title('\rm Inversion with r(w)=1');
t2pos = get(t2,'position');
set(t2,'position',t2pos+[0 0.01 0]);
set(gca,'box','off')
ylim([0.5 1])

% load abs results

load(['figures_paper/lambda_abs_',low,'_',averaging,'.mat'])

%load(['figures_paper/lambda_abs_',low,'_',averaging,'_rms1em5.mat'])


surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808, ...
310.7297, 316.1383];

subplot(2,2,3)

plot(surf_temp,lambda_turb_omega,'b','linewidth',1.2); hold on;
plot(surf_temp,lambda_turb_omega_QG,'b--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega,'r','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega_QG,'r--');
x1 = xlabel('Global-mean surface air temperature (K)');
x1pos = get(x1,'position');
set(x1,'position',x1pos+[0 -0.03 0]);
y2 = ylabel('Asymmetry parameter $\lambda$');
%y2 = ylabel(sprintf('Asymmetry parameter \\lambda %g\t'));
y2pos = get(y2,'position');
set(y2,'position',y2pos+[-5 0 0]);
set(gca,'FontSize',11);
%title('\rm RHS+r') 
ylim([0.5 1])
set(gca,'box','off')

subplot(2,2,4)

plot(surf_temp,lambda_turb_omega,'b','linewidth',1.2); hold on;
plot(surf_temp,lambda_turb_omega_QG_rhs,'b--','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega,'r','linewidth',1.2); hold on;
plot(surf_temp,lambda_mode_omega_QG_rhs,'r--','linewidth',1.2);
x2 = xlabel('Global-mean surface air temperature (K)');
x2pos = get(x2,'position');
set(x2,'position',x2pos+[0 -0.03 0]);
set(gca,'box','off')
set(gca,'box','off')
%ylabel('\lambda');
set(gca,'FontSize',11);
%title('\rm RHS') 
ylim([0.5 1])
set(gca,'box','off')

% % a,b,c,d
% annotation('textbox', [0.07, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.5, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.07, 0.47, 0.01, 0.03], 'String', "(c)",'FontSize',12,'EdgeColor','none')
% annotation('textbox', [0.5, 0.47, 0.01, 0.03], 'String', "(d)",'FontSize',12,'EdgeColor','none')

% a,b,c,d
annotation('textbox', [0.07, 0.94, 0.01, 0.03], 'String', "(a)",'FontSize',12,'EdgeColor','none','Interpreter','latex')
annotation('textbox', [0.5, 0.94, 0.01, 0.03], 'String', "(b)",'FontSize',12,'EdgeColor','none','Interpreter','latex')
annotation('textbox', [0.07, 0.47, 0.01, 0.03], 'String', "(c)",'FontSize',12,'EdgeColor','none','Interpreter','latex')
annotation('textbox', [0.5, 0.47, 0.01, 0.03], 'String', "(d)",'FontSize',12,'EdgeColor','none','Interpreter','latex')

% text boxes

% %annotation('textbox', [0.1, 0.5, 0.03, 0.3], 'String', "Asymmetry parameter \lambda",'EdgeColor','none')
% %annotation('textbox', [0.1, 0.3, 0.03, 0.3], 'String', "Asymmetry parameter \lambda",'EdgeColor','none')
% xlim=get(gca,'XLim');
% ylim=get(gca,'YLim');
% ht = text(0.9*xlim(1)+0.1*xlim(2),0.9*ylim(1)+0.1*ylim(2),'Asymmetry parameter \lambda');
% set(ht,'Rotation',90)
% set(ht,'FontSize',12)

annotation('textbox', [0.91, 0.8, 0.01, 0.03], 'String', "Reduced Stability Simulations",'EdgeColor','none','Interpreter','latex')
annotation('textbox', [0.91, 0.3, 0.01, 0.03], 'String', "Global Warming Simulations",'EdgeColor','none','Interpreter','latex')

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex')

saveas(gcf,['figures_paper/lambda_qg_and_gcm_',low,'_',averaging],'epsc');