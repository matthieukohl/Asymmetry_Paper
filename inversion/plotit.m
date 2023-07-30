figure(1)
plot(surf_temp,Skew_total); hold on;
plot(surf_temp,Skew_Deformation)
xlabel('T_{g}')
ylabel('Skew')
title('Skew J(\tau,del \phi) vs. Deformation at 500hPA for abs-turb.')
legend('Vorticity Advection','Deformation','location','northwest')
legend boxoff 
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_turbulent_1.fig')

figure(2)
%set(gcf,'DefaultAxesFontSize', 14)
plot(surf_temp,Skew_total); hold on;
plot(surf_temp,Skew_shear_u); hold on;
plot(surf_temp,Skew_shear_v); hold on;
plot(surf_temp,Skew_vor); hold on;
%plot(surf_temp,Skew_Product)
xlabel('T_{g}')
ylabel('Skew')
title('Skew J(\tau,del \phi) Decomposition at 500hPA for abs-turb.')
legend('Total','du/dp','dv/dp','Zeta','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_turbulent_2.fig')

figure(3)
plot(surf_temp,Skew_total); hold on;
plot(surf_temp,Skew_A1_bb); hold on;
plot(surf_temp,Skew_A1_pp); hold on;
plot(surf_temp,Skew_A1_pb); hold on;
plot(surf_temp,Skew_A1_bp)
xlabel('T_{g}')
ylabel('Skew')
title('Skew J(\tau,del \phi) Decomposition at 500hPA for abs-turb.')
legend('A1','A1_{bb}','A1_{pp}','A1_{pb}','A1_{bp}','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/skew_abs_turbulent_3.fig')

figure(4)
%set(gcf,'DefaultAxesFontSize', 14)
semilogy(surf_temp,Norm_A1); hold on;
semilogy(surf_temp,Norm_A3); hold on;
semilogy(surf_temp,Norm_Deformation); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('RHS-terms at 500hPA for abs-turb.')
legend('Vort. Advection','Beta','Deformation','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_turbulent_1.fig')

figure(5)
%set(gcf,'DefaultAxesFontSize', 12)
%semilogy(surf_temp,Norm_zeta); hold on;
semilogy(surf_temp,Norm_zeta_bar); hold on;
semilogy(surf_temp,Norm_zeta_prime); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('Zeta Decomposition at 500hPA for abs-turb.')
legend('Bar','Prime','location','northwest')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_turbulent_2.fig')

figure(6)
%set(gcf,'DefaultAxesFontSize', 14)
%semilogy(surf_temp,Norm_du_dp); hold on;
semilogy(surf_temp,Norm_du_dp_bar); hold on;
semilogy(surf_temp,Norm_du_dp_prime); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('Shear u-Decomposition at 500hPA for abs-turb.')
legend('Bar','Prime','location','northeast')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_turbulent_3.fig')

figure(7)
%set(gcf,'DefaultAxesFontSize', 14)
%semilogy(surf_temp,Norm_du_dp); hold on;
semilogy(surf_temp,Norm_dv_dp_bar); hold on;
semilogy(surf_temp,Norm_dv_dp_prime); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('Shear v-Decomposition at 500hPA for abs-turb.')
legend('Bar','Prime','location','northeast')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_turbulent_4.fig')

figure(8)
%set(gcf,'DefaultAxesFontSize', 14)
semilogy(surf_temp,Norm_A1); hold on;
semilogy(surf_temp,Norm_A1_bb); hold on;
semilogy(surf_temp,Norm_A1_pp); hold on;
semilogy(surf_temp,Norm_A1_pb); hold on;
semilogy(surf_temp,Norm_A1_bp); hold on;
xlabel('T_{g}')
ylabel('Norm')
title('A Decomposition at 500hPA for abs-turb.')
legend('Total','A1_{bb}','A1_{pp}','A1_{bp}','A1_{pb}','location','northeast')
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
legend boxoff 
set(gca,'linewidth',1.5)
savefig('/net/aimsir/archive1/mkohl/vorticity/diagnostic_figures/mag_abs_turbulent_5.fig')
