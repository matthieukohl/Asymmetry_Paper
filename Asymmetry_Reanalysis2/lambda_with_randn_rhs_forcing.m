% plot lambda with randon noise RHS forcing

R = logspace(-3,0,200);
%alpha = [-1,0,1/4,1,2];
alpha = [-1/2,-4/5];

trials = 10;

asym_new = zeros(length(R),length(alpha));

for ii = 1:length(R)
    for jj = 1:length(alpha)
        for kk = 1:trials
        
       asym = toy_model_randn(R(ii),alpha(jj));
       asym_new(ii,jj) = asym_new(ii,jj)+asym;
       
        end
        
    end
end
    
asym_new = asym_new/trials;

% load lambda from the GCM and Inversions

low = 'w';

averaging = 'zm';

load(['/net/aimsir/archive1/mkohl/inversion/figures_paper/lambda_r_',low,'_',averaging,'.mat'])

r_gcm = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

grayColor = [.7 .7 .7];

figure(1)
semilogx(R,asym_new); hold on;
semilogx(r_gcm,lambda_turb_omega,'Color',grayColor); hold on;
semilogx(r_gcm,lambda_turb_omega_QG,'Color',grayColor);
xlabel('Reduction factor r')
%legend('\alpha=-1','\alpha=0','\alpha=0.25','\alpha=1','alpha=2');
legend('\alpha = -1/2','\alpha=-4/5');
legend boxoff


% test the arbssnoise function (noise with spectral slope alpha)

L = 12*pi;
N = 512;

alpha = -1/4;

k = 2*pi/L*[0:N/2,-N/2+1:-1];

noise_power = 0;

for ii = 1:100

noise = arbssnoise(N,alpha);
noise_power = noise_power + abs(fft(noise)).^2;

end

noise_power = noise_power/100;
noise_power = noise_power/max(noise_power);

figure
loglog(k,noise_power,'b'); hold on;
loglog(k,k.^alpha,'k');
xlabel('wavenumber');
ylabel('power')