% compare the seasonal cycle based on w vs. rhs spectrum

close all;
clear;

load('fields_w_theory.mat','lambda_era5_NH_xy_weighted','lambda_NH_theory','lambda_NH_modal_theory','lambda_era5_SH_xy_weighted','lambda_SH_theory','lambda_SH_modal_theory','r_500_NH','r_500_SH','k_NH_500','LD_NH_500','k_SH_500','LD_SH_500')

lambda_NH_w = lambda_era5_NH_xy_weighted;
lambda_SH_w = lambda_era5_SH_xy_weighted;

lambda_NH_theory_w = lambda_NH_theory;
lambda_SH_theory_w = lambda_SH_theory;

lambda_NH_modal_theory_w = lambda_NH_modal_theory;
lambda_SH_modal_theory_w = lambda_SH_modal_theory;

r_500_NH_w = r_500_NH;
r_500_SH_w = r_500_SH;

k_NH_w = k_NH_500.*LD_NH_500;
k_SH_w = k_SH_500.*LD_SH_500;

clear('lambda_era5_NH_xy_weighted','lambda_era5_SH_xy_weighted','lambda_NH_theory','lambda_SH_theory','r_500_NH','r_500_SH');

load('fields_rhs_theory.mat','lambda_era5_NH_xy_weighted','lambda_NH_theory','lambda_NH_modal_theory','lambda_era5_SH_xy_weighted','lambda_SH_theory','lambda_SH_modal_theory','r_500_NH','r_500_SH','k_NH_rhs','LD_NH_500','k_SH_rhs','LD_SH_500')

lambda_NH_rhs = lambda_era5_NH_xy_weighted;
lambda_SH_rhs = lambda_era5_SH_xy_weighted;

lambda_NH_theory_rhs = lambda_NH_theory;
lambda_SH_theory_rhs = lambda_SH_theory;

lambda_NH_modal_theory_rhs = lambda_NH_modal_theory;
lambda_SH_modal_theory_rhs = lambda_SH_modal_theory;

r_500_NH_rhs = r_500_NH;
r_500_SH_rhs = r_500_SH;

k_NH_rhs = k_NH_rhs.*LD_NH_500;
k_SH_rhs = k_SH_rhs.*LD_SH_500;

% make figure

set(groot,'DefaultTextInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])
%x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)

h1 = plot(lambda_era5_NH_xy_weighted,'linewidth',1.5,'Color','b'); hold on;
%h11 = plot(lambda_era5_NH_x,'linewidth',1.5,'Color','b','Linestyle','--'); hold on;

h2 = plot(lambda_NH_theory_w,'linewidth',1.5,'Color','k','Linestyle','-.'); hold on;
h3 = plot(lambda_NH_theory_rhs,'linewidth',1.5,'Color','k','Linestyle',':'); hold on;

h4 = plot(lambda_NH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--');hold on
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
l1 = legend([h4,h3,h2,h1],{'1-D Modal Theory','1-D Toy Model rhs','1-D Toy Model w','ERA5 Coarse-Grained'},'Location','NorthWest'); legend boxoff
l1pos = get(l1,'position');
set(l1,'position',l1pos+[-0.01 0.02 0 0]);
l1.FontSize = 9;
title(['$ \lambda$ (',num2str(500),'hPa): ',num2str(30),'$^{\circ}$-',num2str(70),'$^{\circ}$',' NH'])
set(gca,'box','off')

subplot(1,2,2)

h1 = plot(lambda_era5_SH_xy_weighted,'linewidth',1.5,'Color','b'); hold on;
%h11 = plot(lambda_era5_SH_x_weighted,'linewidth',1.5,'Color','b','Linestyle','--'); hold on;
h2 = plot(lambda_SH_theory_w,'linewidth',1.5,'Color','k','Linestyle','-.'); hold on;
h3 = plot(lambda_SH_theory_rhs,'linewidth',1.5,'Color','k','Linestyle',':'); hold on;
h4 = plot(lambda_SH_modal_theory,'linewidth',1.5,'Color','k','Linestyle','--'); hold on;
ylabel('Asymmetry parameter $\lambda$'); ylim([0.5 0.85])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
l2 = legend([h4,h3,h2,h1],{'1-D Modal Theory','1-D Toy Model rhs','1-D Toy Model w','ERA5 Coarse-Grained'},'Location','NorthWest'); legend boxoff
l2pos = get(l2,'position');
set(l2,'position',l2pos+[-0.01 0.02 0 0]);
l2.FontSize = 9;
title(['$ \lambda$ (',num2str(500),'hPa): ',num2str(30),'$^{\circ}$-',num2str(70),'$^{\circ}$',' SH'])
set(gca,'box','off')

annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

%saveas(gcf,'/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_seasonal_cycle_reanalysis_vs_theory_xy_mean','epsc');


figure(2);

plot(k_NH_w,'k','Linewidth',1.5); hold on;
plot(k_NH_rhs,'k--','Linewidth',1.5); hold on;
plot(k_SH_w,'b','Linewidth',1.5); hold on;
plot(k_SH_rhs,'b--','Linewidth',1.5)
legend('NH w','NH rhs','SH w','SH rhs','Location','NorthWest'); legend boxoff;
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',12)
set(gca,'box','off')
ylabel('Wavenumber k')
%yticks([0.1 0.2 0.3 0.4 0.5])
title('Wavenumber k (500hPa)')

