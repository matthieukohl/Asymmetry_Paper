% compare lambda(r) 1d vs 2d for different k

close all; clear;

r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];
k = [1.2,1.8];

lambda_1d_vs_r = zeros(length(r),length(k));
lambda_2d_vs_r = zeros(length(r),length(k));




for ii = 1:length(r)
    for jj = 1:length(k)
L = 2*pi/k(jj);
Nx = 300;
Ny = 300;
x = linspace(0,L,Nx); 
y = linspace(0,L,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
[X,Y] = meshgrid(x,y);
RHS = sin(k(jj).*X).*sin(k(jj)*Y);
     
% Inversions  
%[~,lambda_1d_vs_r(ii,jj)] = Omega_Solver1D(RHS(:),r(ii),Nx*Nx,dx);
[~,lambda_1d_vs_r(ii,jj)] = Omega_Solver(RHS,r(ii),Nx,Ny,dx,dy);


% double resolution
L = 2*pi/k(jj);
Nx = 400;
Ny = 400;
x = linspace(0,L,Nx); 
y = linspace(0,L,Ny);
dx = x(2)-x(1);
dy = y(2)-y(1);
[X,Y] = meshgrid(x,y);
RHS = sin(k(jj).*X).*sin(k(jj)*Y);

[~,lambda_2d_vs_r(ii,jj)] = Omega_Solver(RHS,r(ii),Nx,Ny,dx,dy);


    end
end

% load GCM fields

% choose lower boundary condition

%low = 'w0';
low = 'w';

% choose averaging technique: 'zm' zonal mean, 'wzm' without zonal mean

averaging = 'zm';

% % change interpreter to latex
% 
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% load r simulation results

load(['figures_paper/lambda_r_',low,'_',averaging,'.mat'])

r_gcm = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];


figure
%semilogx(r,lambda_1d_vs_r(:,1),'k','Linewidth',1.6); hold on;
%semilogx(r,lambda_2d_vs_r(:,1),'r','Linewidth',1.6); hold on;
semilogx(r,lambda_1d_vs_r(:,2),'r','Linewidth',1.6); hold on;
semilogx(r,lambda_2d_vs_r(:,2),'r--','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega,'b','Linewidth',1.6); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'b--','Linewidth',1.6);
xlabel('Reduction factor r');
ylabel('Asymmetry parameter \lambda');