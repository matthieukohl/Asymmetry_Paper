% plot the toy-model results

close all; clear;

% define r

r = linspace(1e-3,1,300);
%r = logspace(-7,0,100);

% define k 

%k = 1:0.3:1.9;
%k = 1.3:0.3:2.2;

%k = 1.2:0.3:2.1;
k = [1.2,1.5,1.8];
k = [k,6.6];

%k = [1.8,2.3,6.6,11];


% asymmetry

lambda = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
     [~,lambda(ii,jj)] = toy_model(r(ii),k(jj));
        
    end
end

% define gradient colormap
%grad=colorGradient([0.85 0.325 0.098],[0.9290 0.6940 0.1250],size(lambda,2));
%grad=colorGradient([1 0.7812 0.4975],[0 0 0],size(lambda,2));
grad = copper(size(lambda,2));
%grad=colorGradient([0 0 1],[1 0 0],size(lambda,2));

grad = flip(grad,1);

figure(1)
x0=10; y0=10; width=900; height=300; set(gcf,'position',[x0,y0,width,height])

%x0=10; y0=10; width=400; height=300; set(gcf,'position',[x0,y0,width,height])

subplot(1,2,1)
h = semilogx(r,lambda,'Linewidth',1.6); hold on;
set(h,{'color'},num2cell(grad,2)); % apply gradient colormap
%legend('k=1.0','k=1.3','k=1.6','k=1.9','k=6.0');
%legend('k=1.2','k=1.5','k=1.8','k=2.1','k=6.0');
legend('k=1.2','k=1.5','k=1.8','k=6.6');
%legend('k=1.8','k=2.3','k=6.6','k=11');
%legend('k=1.2','k=1.5','k=1.8');
legend boxoff
set(gca,'FontSize',12)
xlabel('Reduction factor r')
ylabel('Asymmetry parameter $\lambda$')
%title('Asymmetry parameter $\lambda$')
ylim([0.5 1])
set(gca,'box','off')

%saveas(gcf,'figures_paper/toy_model_asymmetry_and_w_profile_alternative_one_panel','epsc');

% calculate the vertical velocity profile at fixed r

%[~,k1p6] = min(abs(k-1.6));
[~,k1p6] = min(abs(k-1.8));

r_sample = [0.01,0.1,0.5,1.0];

w_profile = zeros(300,length(r_sample));

for ii = 1:length(r_sample)
[w_profile(:,ii),~] = toy_model(r_sample(ii),k(k1p6));
end

x = linspace(0,2*pi/k(k1p6),300);

% define gradient colormap
%grad=colorGradient([0.85 0.325 0.098],[0.9290 0.6940 0.1250],size(w_profile,2));
grad = copper(size(w_profile,2));
grad = flip(grad,1);

subplot(1,2,2)
h = plot(x,w_profile,'Linewidth',1.6); 
set(h,{'color'},num2cell(grad,2)); % apply gradient colormap
legend('r=0.01','r=0.1','r=0.5','r=1.0','Location','NorthWest');
legend boxoff
set(gca,'FontSize',12)
xlabel('x')
%ylabel('Vertical velocity w')
%title('w profile (k=1.8)')
set(gca,'box','off')
xlim([x(1) x(end)])

%annotation('textbox', [0.06, 0.94, 0.01, 0.03], 'String', "(a)",'EdgeColor','none')
%annotation('textbox', [0.50, 0.94, 0.01, 0.03], 'String', "(b)",'EdgeColor','none')

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%saveas(gcf,'figures_paper/toy_model_asymmetry_and_w_profile_alternative','epsc');
%saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/inversion_abs_rhs','epsc');


% plot lambda vs k for fixed r

% asymmetry
r = [0.001,0.01,0.1,0.4];
k = linspace(1,11,20);

lambda = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
     [~,lambda(ii,jj)] = toy_model(r(ii),k(jj));
        
    end
end

[a,ind] = min(abs(r-0.01));

plot(k,lambda,'Linewidth',1.4);
xlabel('Wavenumber k')
ylabel('Asymmetry parameter $\lambda$');
lgd = legend('r=0.001','r=0.01','r=0.1','r=0.4','Location','NorthWest'); 
legend boxoff;
lgd.FontSize = 11;
xlim([1 11]);
ylim([0.5 1])
set(gca,'box','off')
set(gca,'FontSize',12);

%saveas(gcf,'figures_paper/toy_model_asymmetry_vs_k','epsc');