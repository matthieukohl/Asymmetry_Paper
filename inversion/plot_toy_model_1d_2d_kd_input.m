% plot the 1d and 2d sine and box model where
% we use kd from the centroid of the spec(rhs)/(1+k^2)^2

close all; clear;

r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

kd = [0.48];
kd2 = kd.^2;

% 2d theory

[Kd2,R] = meshgrid(kd2,r);

Ku2 = (-1+sqrt(1+8*R.*Kd2.*(1+2*Kd2)))./(4*R); % 1st root

lambda_2d = 1./(1+Kd2./Ku2);

% 1d theory

lambda_1d = zeros(size(lambda_2d));

for ii = 1:length(r)
    for jj = 1:length(kd)

x0 = 1;
ku = fsolve(@(x)root_1d_toy_model(x,kd(jj),r(ii)),x0);

lambda_1d(ii,jj) = 1/(1+kd/ku);
    
    end

end


lambda_1d_sin = zeros(size(lambda_1d));
lambda_2d_sin = zeros(size(lambda_1d));



for ii = 1:length(kd)
    for jj = 1:length(r)
        
% calculate k for 1d

x0 = 1;
ku = fsolve(@(x)root_1d_toy_model(x,kd(ii),r(jj)),x0);

k = 2*ku*kd/(ku+kd);

[a,lambda_1d_sin(jj,ii)] = toy_model(r(jj),k);

% calculate k for 2d

kd2 = kd^2;
ku2 = (-1+sqrt(1+8*r(jj).*kd2.*(1+2*kd2)))./(4*r(jj)); % 1st root
ku = sqrt(ku2);

k = 2*ku*kd/(ku+kd);

[a,lambda_2d_sin(jj,ii)] = toy_model_2d_alternative(r(jj),k);

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


set(groot,'DefaultTextInterpreter','default');
set(groot,'DefaultLegendInterpreter','default');

figure
semilogx(r,lambda_1d,'b-o'); hold on;
semilogx(r,lambda_1d_sin,'c-o'); hold on;
semilogx(r,lambda_2d,'g-o'); hold on;
semilogx(r,lambda_2d_sin,'m-o');
semilogx(r_gcm,lambda_turb_omega,'k'); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'k--')
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('1d','1d sine','2d','2d sine','GCM','3D Inversion','Location','SouthWest'); legend boxoff;
%legend('1d','2d','2d refined','1d sine','2d sine','2d sin alt','2d sine wp','GCM','3D Inversion'); legend boxoff;
%legend('1d','2d','2d+asym','1d sine','2d sine','GCM','3D Inversion'); legend boxoff;
set(gca,'FontSize',12)
title(['\rm Wavenumber k_d =',num2str(kd)]);
ylim([0.5 0.9])