% plot lambda for w and w_qg for r and abs simulations

clear; close all;

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

figure(1)
%x0=10; y0=10; width=900; height=600; set(gcf,'position',[x0,y0,width,height])

%subplot(2,2,1)
semilogx(r,lambda_turb_omega,'b','linewidth',1.4); hold on;
semilogx(r,lambda_turb_omega_QG,'b--','linewidth',1.4); hold on;
semilogx(r,lambda_mode_omega,'r','linewidth',1.4); hold on;
semilogx(r,lambda_mode_omega_QG,'r--','linewidth',1.4); hold on;
x1 = xlabel('Reduction factor r');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
y1 = ylabel('Asymmetry parameter $\lambda$');
L(1) = plot(nan, nan, 'k');
L(2) = plot(nan, nan, 'k--');
legend(L, {'Simulation', 'Moist Omega Equation'},'box','off','FontSize',12)
lgd.TextColor = 'black';

set(gca,'box','off')
ylim([0.5 1])
%set(y1,'position',y1pos+[-0.0023 0.05 0]);
set(gca,'FontSize',12);


saveas(gcf,['figures_committee_meeting_aug26th/lambda_qg_vs_gcm_r'],'epsc');

figure(2)
semilogx(r,lambda_turb_omega,'b','linewidth',1.4); hold on;
semilogx(r,lambda_turb_omega_QG_rhs,'b--','linewidth',1.4); hold on
semilogx(r,lambda_mode_omega,'r','linewidth',1.4); hold on;
semilogx(r,lambda_mode_omega_QG_rhs,'r--','linewidth',1.4);
y1 = ylabel('Asymmetry parameter $\lambda$');
x1 = xlabel('Reduction factor r');
%x1pos = get(x1,'position');
%set(x1,'position',x1pos+[0 -0.03 0]);
%ylabel('\lambda')
set(gca,'FontSize',12);
L(1) = plot(nan, nan, 'k');
L(2) = plot(nan, nan, 'k--');
legend(L, {'Simulation', 'Moist Omega Equation with r=1'},'box','off','FontSize',12)
lgd.TextColor = 'black';
%lgd.FontSize = 20;
set(gca,'box','off')
ylim([0.5 1])

saveas(gcf,['figures_committee_meeting_aug26th/lambda_qg_vs_gcm_r_rhs'],'epsc');