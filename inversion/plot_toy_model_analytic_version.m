% plot the toy model result based on the semi-analytic method
clear; close all; 

r = 0.01; 
k = 1;

b = fsolve(@(x)root_finder(k,r,x),1);

cd = -sin(k*b)/(1+k^2)*1/cosh(b+pi/(2*k));
cu = -sin(k*b)/(1+r*k^2)*1/cosh(1/sqrt(r)*(b-pi/(2*k)));

wd = @(x) -sin(k*x)/(1+k^2)+cd*cosh(x-pi/(2*k));
wu = @(x) -sin(k*x)/(1+r*k^2)+cu*cosh(1/sqrt(r)*(x-3*pi/(2*k)));

xd = linspace(-b,pi/k+b,500);
xu = linspace(pi/k+b,2*pi/k-b,500);

xd2_plot = linspace(2*pi/k-b,2*pi/k,100);
xd2_cal = linspace(-b,0,100);
x = [xd,xu];

figure(1)
plot(xd,wd(xd),'b'); hold on;
plot(xu,wu(xu),'b');
set(gca,'FontSize',12)
xlabel('x')
ylabel('w')
xlim([x(1) x(end)])
set(gca,'box','off');

% check mass conservation

Int = integral(wd,xd(1),xd(end))+integral(wu,xu(1),xu(end));

disp('int=')
Int

% compare to numerical inversions

[w_inversion,~] = toy_model(r,k);
x = linspace(0,2*pi/k,300);

% compare asymmetries
disp('\lambda inversion')
asymmetry(w_inversion)
disp('\lambda numerical')
asymmetry_integral_version([xd,xu],[wd(xd),wu(xu)])

figure(2)
plot(x,w_inversion,'b'); hold on;
plot(xd,wd(xd),'k--','Linewidth',1.4); hold on;
plot(xu,wu(xu),'k--','Linewidth',1.4); hold on;
plot(xd2_plot,wd(xd2_cal),'k--','Linewidth',1.4);
legend('w_{inv}','w_{analytic}'); 
legend boxoff
set(gca,'FontSize',12)
xlabel('x')
ylabel('w')
xlim([x(1) x(end)])
set(gca,'box','off');

% make plots of lambda(r,k)

r = logspace(-3,0,100);

% define k 

k = 1:0.3:1.9;
k = [k,6];

% asymmetry

lambda = zeros(length(r),length(k));
lambda_numerical = zeros(length(r),length(k));
bb = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
rr = r(ii); kk = k(jj);

b = fsolve(@(x)root_finder(kk,rr,x),0.1);

bb(ii,jj) = b;

cd = -sin(kk*b)/(1+kk^2)*1/cosh(b+pi/(2*kk));
cu = -sin(kk*b)/(1+rr*kk^2)*1/cosh(1/sqrt(rr)*(b-pi/(2*kk)));

wd = @(x) -sin(kk*x)/(1+kk^2)+cd*cosh(x-pi/(2*kk));
wu = @(x) -sin(kk*x)/(1+rr*kk^2)+cu*cosh(1/sqrt(rr)*(x-3*pi/(2*kk)));

xd = linspace(-b,pi/kk+b,300);
xu = linspace(pi/kk+b,2*pi/kk-b,300);

w = [wd(xd),wu(xu)];

%lambda(ii,jj) = asymmetry(w);
lambda(ii,jj) = asymmetry_integral_version([xd,xu],[wd(xd),wu(xu)]);

[~,lambda_numerical(ii,jj)] = toy_model(r(ii),k(jj));

if ii==1
    
% test asymptotic solutions

b = pi./(2*kk)-1./kk.^2.*tanh(pi./kk);   
%b = 1./tanh(pi./(2*k)); 

cd = -sin(kk*b)/(1+kk^2)*1/cosh(b+pi/(2*kk));
cu = -sin(kk*b)/(1+rr*kk^2)*1/cosh(1/sqrt(rr)*(b-pi/(2*kk)));

wd = @(x) -sin(kk*x)/(1+kk^2)+cd*cosh(x-pi/(2*kk));
wu = @(x) -sin(kk*x)/(1+rr*kk^2)+cu*cosh(1/sqrt(rr)*(x-3*pi/(2*kk)));

xd = linspace(-b,pi/kk+b,300);
xu = linspace(pi/kk+b,2*pi/kk-b,300);

w = [wd(xd),wu(xu)];

%lambda(ii,jj) = asymmetry(w);
lambda_asym(jj) = asymmetry_integral_version([xd,xu],[wd(xd),wu(xu)]);
end
    

    end
end

%b_asym = pi./(2*k)-1./k.^2.*tanh(pi./k);
%lambda_asym = 1-(pi./k-2*b_asym)./(3*pi./k-2*b_asym);

figure(3)
%x0=10; y0=10; width=00; height=300; set(gcf,'position',[x0,y0,width,height])

semilogx(r,lambda,'Linewidth',1.6); hold on;
semilogx(r,lambda_numerical,'k--','Linewidth',1.6);
scatter(r(1)*ones(length(k),1),lambda_asym)
legend('k=1.0','k=1.3','k=1.6','k=1.9','k=6.0');
legend boxoff
set(gca,'FontSize',12)
xlabel('Reduction factor r')
ylabel('Asymmetry parameter $\lambda$')
title('Asymmetry')
ylim([0.5 1])
set(gca,'box','off')

% check accuracy of asymptotic solution

r = 1e-5;

% define k 

k = 1:0.3:1.9;
k = [k,6];

% asymmetry

lambda = zeros(length(r),length(k));
lambda_numerical = zeros(length(r),length(k));
bb = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
rr = r(ii); kk = k(jj);

b = fsolve(@(x)root_finder(kk,rr,x),0.5);

bb(ii,jj) = b;


    end
end

b_asym = pi./(2*k)-1./k.^2.*tanh(pi./k);
%eps = 1./(2*k.^2).*(-sech(pi./k).^2+sqrt(sech(pi./k).^4-4*k.^2.*tanh(pi./k)));
%b_asym = 1./tanh(pi./(2*k));


%b_asym = pi./(2*k)+eps;

plot(k,bb); hold on;
plot(k,b_asym);
