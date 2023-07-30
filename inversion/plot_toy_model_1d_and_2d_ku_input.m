% plot the 1d and 2d sine and box model where
% we use ku from the centroid of the w spectrum as input

close all; clear;

r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

ku = [1.8];
ku2 = ku.^2;

% 2d theory

[Ku2,R] = meshgrid(ku2,r);

Kd2 = 1/4*(-1 + sqrt(1+8*Ku2+16*R.*Ku2.^2));

lambda_2d = 1./(1+Kd2./Ku2);

% 1d theory

lambda_1d = zeros(size(lambda_2d));

for ii = 1:length(r)
    for jj = 1:length(ku)

x0 = 1;
kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku(jj),r(ii)),x0);

lambda_1d(ii,jj) = 1/(1+kd/ku);
    
    end

end


lambda_1d_sin = zeros(size(lambda_1d));
lambda_2d_sin = zeros(size(lambda_1d));



for ii = 1:length(ku)
    for jj = 1:length(r)
        
% calculate k for 1d

x0 = 1;
kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku(ii),r(jj)),x0);

k = 2*ku*kd/(ku+kd);

[a,lambda_1d_sin(jj,ii)] = toy_model(r(jj),k);

% calculate k for 2d

ku2 = ku^2;
kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r(jj)*ku2^2));

kd = sqrt(kd2);

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
title(['\rm Wavenumber k_u =',num2str(kd)]);
ylim([0.5 0.9])