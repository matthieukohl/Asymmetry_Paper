% Compare different theoretical predictions for values taken at 500hPa,mid,
% or average

close all; clear;

load('data_theory/theory_NH_complete2');
load('data_theory/theory_SH_complete2');

hfig = figure(1);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(lambda_NH_theory_mid,'Linewidth',1.5); hold on
plot(lambda_NH_theory_av,'Linewidth',1.5); hold on;
plot(lambda_NH_theory_500,'Linewidth',1.5); hold on;
ylabel('\lambda'); ylim([0.5 0.9]) 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('Mid','Av','500','Location','Northwest'); legend boxoff
title(['Northern Hemisphere'])

subplot(1,2,2)

plot(lambda_SH_theory_mid,'Linewidth',1.5); hold on
plot(lambda_SH_theory_av,'Linewidth',1.5); hold on;
plot(lambda_SH_theory_500,'Linewidth',1.5); hold on;
ylabel('\lambda'); ylim([0.5 0.9])
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
title(['Southern Hemisphere'])

hfig = figure(2);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(r_mid_NH,'Linewidth',1.5); hold on
plot(r_av_NH,'Linewidth',1.5); hold on;
plot(r_500_NH,'Linewidth',1.5); hold on;
ylabel('r'); 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('Mid','Av','500','Location','Northwest'); legend boxoff
title(['Northern Hemisphere'])

subplot(1,2,2)

plot(r_mid_SH,'Linewidth',1.5); hold on
plot(r_av_SH,'Linewidth',1.5); hold on;
plot(r_500_SH,'Linewidth',1.5); hold on;
ylabel('r'); 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
title(['Southern Hemisphere'])

hfig = figure(3);
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 2 1])

subplot(1,2,1)

plot(k_NH_mid.*L_mid_NH,'Linewidth',1.5); hold on
plot(k_NH_av.*L_av_NH','Linewidth',1.5); hold on;
plot(k_NH_500.*L_500_NH,'Linewidth',1.5); hold on;
ylabel('r'); 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('Mid','Av','500','Location','Northwest'); legend boxoff
title(['Northern Hemisphere'])

subplot(1,2,2)

plot(-k_SH_mid.*L_mid_SH,'Linewidth',1.5); hold on
plot(-k_SH_av.*L_av_SH,'Linewidth',1.5); hold on;
plot(-k_SH_500.*L_500_SH,'Linewidth',1.5); hold on;
ylabel('r'); 
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
%legend('NCEP2(x)','NCEP2(xt)','ERAI(x)','ERAI(xt)','ERA5(x)','ERA5(xt)','Turb','Modal','ZG','PG','Location','Northeast'); legend boxoff
title(['Southern Hemisphere'])

