% check if the 2d inversion with RHS = sin(kX)
% instead of RHS = sin(kX)*sin(kY) agrees with
% the 1d inversion RHS = sin(kx)

close all; clear;

N = 300;
k = 1.8;
x = linspace(0,2*pi/k,N);

r = [0.01,0.05,0.1,0.2,0.4,0.6,0.8,1.0];

% calculations
w_1d_vs_r = zeros(N,length(r));
w_2d_vs_r = zeros(N,N,length(r));

lambda_1d = zeros(size(r));
lambda_2d = zeros(size(r));

for ii = 1:length(r)

[w_1d_vs_r(:,ii),lambda_1d(ii)] = toy_model(r(ii),k);
[w_2d_vs_r(:,:,ii),lambda_2d(ii)] = toy_model_2d_alternative_zonal(r(ii),k);

end

% compare the w-profiles at r=0.01
w_1d = w_1d_vs_r(:,1);
w_2d = w_2d_vs_r(:,:,1);

[C,I] = max(w_2d(:));
[I1,I2] = ind2sub(size(w_2d),I);

w_2d_slice = squeeze(w_2d(I1,:));

figure
plot(x,w_1d,'b'); hold on;
plot(x,w_2d_slice,'r--');
xlabel('x')
ylabel('w')
set(gca,'FontSize',12);
legend('sin(kx)','sin(kX)','Location','NorthWest'); legend boxoff;

% compare the asymmetries

figure

semilogx(r,lambda_1d,'b'); hold on;
semilogx(r,lambda_2d,'r--');
xlabel('Reduction factor r');
ylabel('Asymmetry factor \lambda');
legend('sin(kx)','sin(kX)'); legend boxoff
set(gca,'FontSize',12)

