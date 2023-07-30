% test the 2d box model vs the 2d sine model

%r = linspace(0.01,1,10);
r = logspace(-2,0,10);
k = 1.8;

% 1d calculation

lambda_box_1d = zeros(length(r),length(k));
lambda_sine_1d = zeros(length(r),length(k));

for ii = 1:length(r)
    
lambda_box_1d(ii) = toy_model_1d_non_sine(r(ii),k);
lambda_sine_1d(ii) = toy_model(r(ii),k,1);

end

% 2d calculation

lambda_box = zeros(length(r),length(k));
lambda_sine = zeros(length(r),length(k));

for ii = 1:length(r)
    
lambda_box(ii) = toy_model_2d_non_sine(r(ii),k);
lambda_sine(ii) = toy_model_2d_alternative(r(ii),k);

end

figure

semilogx(r,lambda_box_1d,'b'); hold on;
semilogx(r,lambda_sine_1d,'b--'); hold on;
semilogx(r,lambda_box,'r'); hold on;
semilogx(r,lambda_sine,'r--');
xlabel('Reduction factor r')
ylabel('Asymmetry parameter \lambda')
legend('box 1d','sine 1d','box 2d','sine 2d'); legend boxoff;
set(groot,'DefaultTextInterpreter','default');
set(groot,'DefaultLegendInterpreter','default');
set(gca,'FontSize',12)
title(['\rm Wavenumber k=',num2str(k)])
