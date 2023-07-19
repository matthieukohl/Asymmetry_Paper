% plot the toy model:
% 2ku^2*r*wu + wu = - (2kd^2 * wd + wd)
% wu/wd = ku^2 / kd^2
% that arise from moist omega and mass conservation equation

close all; clear;

% kd = [1.8];
% kd2 = kd.^2;
r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

ku = [1.8];

ku2 = ku^2;
kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));

kd = sqrt(kd2);

k = 2*ku*kd/(ku+kd);

kd = k;
kd2 = kd.^2;

% 2d calculation

[Kd2,R] = meshgrid(kd2,r);

Ku2 = (-1+sqrt(1+8*R.*Kd2.*(1+2*Kd2)))./(4*R); % 1st root
%Ku2 = (1-sqrt(1+8*R.*Kd2.*(1-2*Kd2)))./(4*R); % 2nd root

%Ku2(1,:) = Kd2(1,:).*(1-2*Kd2(1,:));

lambda_2d = 1./(1+Kd2./Ku2);

% 1d calculation
lambda_1d = zeros(size(lambda_2d));
lambda_2d_refined = zeros(size(lambda_2d));

lambda_1d_sin = zeros(size(lambda_1d));
lambda_1d_sin_area = zeros(size(lambda_1d));
lambda_2d_sin = zeros(size(lambda_1d));
lambda_2d_sin_alt = zeros(size(lambda_1d));
lambda_2d_sin_alt_area = zeros(size(lambda_1d));
lambda_2d_sin_wp = zeros(size(lambda_1d));

x0 = [20,20,20,4,4,4,4,4];

for ii = 1:length(kd)
    for jj = 1:length(r)
ku = fsolve(@(x)root_1d_toy_model(x,kd(ii),r(jj)),x0(jj));
lambda_1d(jj,ii) = 1/(1+kd(ii)/ku);
[a,lambda_1d_sin(jj,ii),lambda_1d_sin_area(jj,ii)] = toy_model(r(jj),kd(ii));
[a,lambda_2d_sin(jj,ii)] = toy_model_2d(r(jj),kd(ii));
[a,lambda_2d_sin_alt(jj,ii),lambda_2d_sin_alt_area(jj,ii)] = toy_model_2d_alternative(r(jj),kd(ii));
[a,lambda_2d_sin_wp(jj,ii)] = toy_model_2d_alternative_without_periodicity(r(jj),kd(ii));
[lambda_2d_refined(jj,ii)] = toy_model_non_sine_2d_refined(r(jj),kd(ii));

    end
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

load('test500.mat');

set(groot,'DefaultTextInterpreter','default');
set(groot,'DefaultLegendInterpreter','default');

figure
semilogx(r,lambda_1d,'b-o'); hold on;
semilogx(r,lambda_1d_sin,'c-o'); hold on;
semilogx(r,lambda_2d_refined,'r-o'); hold on;
semilogx(r,lambda_2d,'g-o'); hold on;
semilogx(r,lambda_2d_sin,'m-o');
semilogx(r,lambda_2d_sin_alt,'m--');
%semilogx(r,lambda_2d_sin_wp,'m--'); hold on;
semilogx(r_gcm,lambda_turb_omega,'k'); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'k--')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('1d','1d sine','2d refined','2d','2d sine old','2d sine new','GCM','3D Inversion','Location','SouthWest'); legend boxoff;
%legend('1d','2d','2d refined','1d sine','2d sine','2d sin alt','2d sine wp','GCM','3D Inversion'); legend boxoff;
%legend('1d','2d','2d+asym','1d sine','2d sine','GCM','3D Inversion'); legend boxoff;
set(gca,'FontSize',12)
title(['\rm Wavenumber k_d =',num2str(kd)]);
ylim([0.5 0.9])

figure
semilogx(r,imag(lambda_1d),'b-o')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
set(gca,'FontSize',12)

% plot for multiple kd

% plot the toy model:
% 2ku^2*r*wu + wu = - (2kd^2 * wd + wd)
% wu/wd = ku^2 / kd^2
% that arise from moist omega and mass conservation equation

close all; clear;

kd = [0.3,0.6,1,1.5,2,5];
kd2 = kd.^2;
r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

% 2d calculation

[Kd2,R] = meshgrid(kd2,r);

%Ku2 = (1+sqrt(1-8*R.*Kd2.*(1-2*Kd2)))./(4*R); % 1st root
%Ku2 = (1-sqrt(1+8*R.*Kd2.*(1-2*Kd2)))./(4*R); % 2nd root

Ku2 = (-1+sqrt(1+8*R.*Kd2.*(1+2*Kd2)))./(4*R); % alternative definition

%Ku2(1,:) = Kd2(1,:).*(1-2*Kd2(1,:));

lambda_2d = 1./(1+Kd2./Ku2);

figure
semilogx(r,lambda_2d,'Linewidth',1.4); hold on;
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('kd=0.3','kd=0.6','kd=1','kd=1.5','kd=2','kd=5','Location','NorthEast'); legend boxoff;
set(gca,'FontSize',12)
title('\rm 2d model')

figure
semilogx(r,imag(lambda_2d),'b-o')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
set(gca,'FontSize',12)

% make a plot

figure

semilogx(r,lambda_1d_sin,'r'); hold on;
semilogx(r,lambda_1d_sin_area,'r--'); hold on;
semilogx(r,lambda_1d,'r-.'); hold on;

semilogx(r,lambda_2d_sin_alt,'b'); hold on;
semilogx(r,lambda_2d_sin_alt_area,'b--'); hold on;
semilogx(r,lambda_2d_refined,'b-.'); hold on;
semilogx(r,lambda_2d,'b:'); hold on;

semilogx(r_gcm,lambda_turb_omega,'k'); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'k--')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('1d sine','1d sine area estimate','1d toy','2d sine','2d sine area estimate','2d toy 2/k = 1/ku + 1/kd','2d toy k = kd','GCM','3D Inversion','Location','SouthWest'); legend boxoff;
%legend('1d','2d','2d refined','1d sine','2d sine','2d sin alt','2d sine wp','GCM','3D Inversion'); legend boxoff;
%legend('1d','2d','2d+asym','1d sine','2d sine','GCM','3D Inversion'); legend boxoff;
set(gca,'FontSize',12)
title(['\rm Wavenumber k_d =',num2str(kd)]);

