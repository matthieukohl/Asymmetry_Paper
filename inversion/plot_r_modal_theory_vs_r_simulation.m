% plot lambda for r modal theory vs r simulation

clear; close all;

% choose lower boundary condition for simulation

%low = 'w0';
low = 'w';

% load r simulation results
load(['figures_paper/lambda_r_',low,'.mat'])

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

% choose boundary condition for theory runs

IC = 'rand';
%IC = 'sinusoidal';

load(['figures_paper/modal_theory_N200_L8pi_',IC,'.mat'])

figure(1)
x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

x = linspace(0,L,N+1);

[a,ind] = min(abs(R-0.1));

subplot(1,2,1)
plot(x,w_final(:,ind),'b','Linewidth',1.8);
xlabel('x');
ylabel('Vertical velocity w');
xlim([x(1) x(end)])
title('\rm w profile (r = 0.1)')
set(gca,'box','off')
set(gca,'FontSize',12);

subplot(1,2,2)
semilogx(r,lambda_mode_omega,'r','linewidth',1.8); hold on;
semilogx(R,lambda,'b','linewidth',1.8); hold on;
x1 = xlabel('Reduction factor r');
y1 = ylabel('Asymmetry parameter $\lambda$');
title('\rm Asymmetry parameter $\lambda$')
legend({'GCM','1-D Modal Theory'},'Location','SouthWest'); legend boxoff
set(gca,'box','off')
ylim([0.5 1])
set(gca,'FontSize',12);



annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

saveas(gcf,'figures_paper/modal_asymmetry_and_w_profile','epsc');