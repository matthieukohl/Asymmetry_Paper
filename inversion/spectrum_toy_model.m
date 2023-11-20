% spectrum toy model: peak at k RHS 
% because L = 2*pi/k, k1 = k_RHS

close all; clear;

k = 1.8;
%L = 2*pi/k;
L = 20;

for ii = 1:1
    
N = 4000*ii;
r = 0.01;


indx = N/2;

x = linspace(0,L,N);
RHS = sin(k*x);

kk = 2*pi/L*[0:N/2,-N/2+1:-1];

[w,lambda] = toy_model_N_L(r,k,N,L);

power = abs(fft(w)).^2;
power = power/max(power);

figure(1)
if ii == 1
semilogx(kk(1:indx),power(1:indx),'b'); hold on;
y1 = get(gca,'ylim');
plot([1.8 1.8],y1,'r');
xlabel('Wavenumber k');
ylabel('power')
set(gca,'FontSize',12);
legend('Spectrum w','Wavenumber k-rhs'); legend boxoff;

loglog(kk(1:indx),power(1:indx),'b'); hold on;
y1 = get(gca,'ylim');
plot([1.8 1.8],y1,'r');
xlabel('Wavenumber k');
ylabel('power')
set(gca,'FontSize',12);
legend('Spectrum w','Wavenumber k-rhs'); legend boxoff;
else 
semilogx(kk(1:indx),power(1:indx),'r--'); hold on;
xlabel('Wavenumber k');
ylabel('power')
set(gca,'FontSize',12)
end

% figure(2)
% plot(x,RHS); hold on;
% plot(x,w)

end


%loglog(kk(1:indx),0.5*kk(1:indx).^-3);

