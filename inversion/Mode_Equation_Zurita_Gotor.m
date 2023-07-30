%% Solves the Non-Dimensional
%% Mode Equation from Zurita-Gotor (2005): rw_xxxx+(1-r*g^2)*w_xx+g^2*w = 0

clear;
close all;

% Define Reduction Factor

R = [1.0;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.18;0.16;0.14;0.12;...
    0.11;0.10;0.08;0.06;0.04;0.02;0.01]; % reduction factor r

lambda = zeros(length(R),1); % initialize lambda factor
sigma = zeros(length(R),1); % initialize growth rate

% Define Grid

N = 200; % number of intervals
L = 4*pi; %4*pi; %2*pi/(sqrt(sqrt(2)-1)); %5*pi % domain length L_max dry unstable mode
x = linspace(0,L,N+1);
dx = L/(N);

% Define Coefficients
% 
% a = 6/dx^4-16/dx^2;
% a1 = 8/dx^2-4/dx^4;
% a2 = 1/dx^4;
% 
% b = -2/(dx^2)-8;
% b1 = 1/(dx^2);

a = 6/dx^4-8/dx^2;
a1 = 4/dx^2-4/dx^4;
a2 = 1/dx^4;

b = -2/(dx^2)-4;
b1 = 1/(dx^2);


% Define the A-matrix 
    
dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1).*ones(1,N-1),-1); % i-1 entries
dAp2 = diag((a2).*ones(1,N-2),2); % i+2 entries
dAm2 = diag((a2).*ones(1,N-2),-2); % i-2 entries

A = (dA + dAp1 + dAm1+ dAp2 + dAm2);

A(1,end) = a1; %periodic BC
A(1,end-1) = a2; %periodic BC
A(2,end) = a2; %periodic BC

A(end,1) = a1; %periodic BC
A(end,2) = a2; %periodic BC
A(end-1,1) = a2; % periodic BC

% Define the B-matrix

dB = diag((b).*ones(1,N)); % diagonal entries
dBp1 = diag((b1).*ones(1,N-1),1); % i+1 entries
dBm1 = diag((b1).*ones(1,N-1),-1); % i-1 entries

B = (dB + dBp1 + dBm1);

B(1,end) = b1; %periodic BC

B(end,1) = b1; %periodic BC

[V,D] = eig(A,B);
[V,D] = sortem(V,D);

sqrt(D(1,1))

w_old = zeros(N,1);
w_new = V(:,1);

% counter 

count = 0;

% Velocity Field

w_field = zeros(N+1,length(R));

for tt = 1:20
    
tt    

count = 0;

%w_new = sin(x(1:end-1)')+sin(7*x(1:end-1)')+sin(5*x(1:end-1)')+sin(10*x(1:end-1)')+sin(3*x(1:end-1)');
    
% Start Iterative Inversion

r = zeros(N,1);

while rms(w_new-w_old)>= 1*10^(-4)
    
count = count + 1;

% Define r coefficient

for jj=1:N
    if w_new(jj)>=0
        r(jj) = R(tt);
    else 
        r(jj) = 1;        
    end
end

% Define Coefficients

a = 6*r/dx^4-8/dx^2;
a1 = 4/dx^2-4*r/dx^4;
a2 = r/dx^4;

b = -2*r/(dx^2)-4;
b1 = 1*r/(dx^2);

% Define the A-matrix 
    
dA = diag((a').*ones(1,N)); % diagonal entries
dAp1 = diag((a1(2:end)').*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1(1:end-1)').*ones(1,N-1),-1); % i-1 entries
dAp2 = diag((a2(3:end)').*ones(1,N-2),2); % i+2 entries
dAm2 = diag((a2(1:end-2)').*ones(1,N-2),-2); % i-2 entries

A = (dA + dAp1 + dAm1+ dAp2 + dAm2);

A(1,end) = a1(end); %periodic BC
A(1,end-1) = a2(end-1); %periodic BC
A(2,end) = a2(end); %periodic BC

A(end,1) = a1(1); %periodic BC
A(end,2) = a2(2); %periodic BC
A(end-1,1) = a2(1); % periodic BC

% Define the B-matrix

dB = diag((b').*ones(1,N)); % diagonal entries
dBp1 = diag((b1(2:end)').*ones(1,N-1),1); % i+1 entries
dBm1 = diag((b1(1:end-1)').*ones(1,N-1),-1); % i-1 entries

B = (dB + dBp1 + dBm1);

B(1,end) = b1(end); %periodic BC

B(end,1) = b1(1); %periodic BC

% Solve eigenvalue problem A*w = sigma^2 B*w

[V,D] = eig(A,B);
[V,D] = sortem(V,D);

diff1 = rms(V(:,1)-w_new);
diff2 = rms(V(:,2)-w_new); 

% w_old = w_new;
% w_new = 0.44*V(:,1)+0.56*w_new;

if diff2>=diff1
w_old = w_new;
w_new = 0.01*V(:,1)+0.99*w_new; % smooth((0.44*V(:,1)+0.56*w_new));

else 
w_old = w_new;
w_new = 0.01*V(:,2)+0.99*w_new;
end

figure(1)
plot(x,[w_new;w_new(1)]); hold on
plot(x,[r;r(1)]);
pause(0.1)
hold off

% xlabel('x');
% ylabel('w');
% title('Most unstable Mode')

% w_1 = V(:,1);
% w_2 = V(:,2);
% w_1 = [w_1;w_1(1)];
% w_2 = [w_2;w_2(1)];
% 
% figure(3)
% plot(x,w_1); hold on;
% plot(x,w_2)
% pause(0.8);
% hold off;
rms(w_new-w_old)
end

w_old = zeros(N,1);

w_final = [w_new;w_new(1)];

w_field(:,tt) = w_final;

lambda(tt) = asymmetry(w_final);

D(1,1)

sigma(tt) = sqrt(D(1,1));


end

% Plot the lambda-profile
% 
theory_lambda = 1./(1+sqrt(R));

lambda_theory_mode = [0.5002, 0.5094, 0.5271, 0.5452, 0.5709, 0.6003, 0.6372, 0.6810,...
0.7448, 0.7547, 0.7733, 0.7882, 0.8063, 0.8169, 0.8296, 0.8520, 0.8755, 0.9007,0.9317, 0.9530];

lambda_theory_correct_N100 = [0.5035,0.5176,0.5313,0.5557,0.5832,0.6154,...
0.6599,0.7067,0.7810,0.7944,0.8079,0.8246,0.8437,0.8476,0.8627,0.8808,...
0.8981,0.9131,0.9513,0.9751];

lambda_theory_correct_N200 = [0.5018,0.5171,0.5352,0.5504,0.5825,0.6114,...
0.6502,0.7085,0.7790,0.7883,0.8112,0.8265,0.8433,0.8512,0.8597,0.8762,...
0.8956,0.9167,0.9415,0.9590];

lambda_theory_correct_8_N100 = [0.4965,0.4968,0.5366,0.5369,0.5467,0.5729,...
0.6383,0.7351,0.7955,0.8246,0.8248,0.8408,0.8585,0.8658,0.8713,0.8906,...
0.9108,0.9261,0.9509,0.9732];

% figure(1)
% plot(x,w_field(:,11),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.16')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(5)
% plot(x,w_field(:,9),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.2')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(6)
% plot(x,w_field(:,8),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.3')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(7)
% plot(x,w_field(:,10),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.18')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(8)
% plot(x,w_field(:,end-2),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.04')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(9)
% plot(x,w_field(:,end-1),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.02')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)
% 
% figure(10)
% plot(x,w_field(:,end-1),'linewidth',1.6); hold on
% xlabel('x','interpreter','latex');
% ylabel('w','interpreter','latex')
% ax = gca;
% ax.XTick = [0,pi,2*pi,3*pi,4*pi];
% ax.XTickLabel = {'0','\pi','2\pi','3\pi','4\pi'};
% title('w-profile for r=0.01')
% set(gca,'fontsize', 14);
% %xlim([0 4*pi])
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'linewidth',1.5)









%saveas(gcf,'Plots_AOFD/w_profile_ZuritaGotor1','epsc')

% figure(2)
% semilogx(R,lambda); hold on;
% semilogx(R,theory_lambda); hold on;
% semilogx(R,lambda_theory_correct_8_N100); hold on;
% xlabel('Updraught Stability Reduction Factor r')
% ylabel('Asymmetry Parameter {\lambda}')
% title('Asymmetry Toy-Model')
% legend('Numerical','Theory','Theory_previous')
%savefig('Zurita_Gotor_Mode1.fig')

%% Make Plots for Generals Mode Equation

 r_factor = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];
% 
% R = [1.0;0.9;0.8;0.7;0.6;0.5;0.4;0.3;0.2;0.18;0.16;0.14;0.12;...
%     0.11;0.10;0.08;0.06;0.04;0.02;0.01];
% 
lambda_omega_mode = [0.9283,0.9173,0.8996,0.8313,0.7820,0.7343,0.6871,...
 0.5674,0.5187,0.4994];

lambda_omegaQG_mode = [0.9193,0.9053,0.8903,0.8200,0.7674,0.7277,0.6756,...
 0.5706,0.5192,0.5043];
% 
% lambda_theory_mode = [0.5002, 0.5094, 0.5271, 0.5452, 0.5709, 0.6003, 0.6372, 0.6810,...
% 0.7448, 0.7547, 0.7733, 0.7882, 0.8063, 0.8169, 0.8296, 0.8520, 0.8755, 0.9007,0.9317, 0.9530];
% 
% 
figure(2)
semilogx(r_factor,lambda_omega_mode,'Color','red'); hold on;
semilogx(r_factor,lambda_omegaQG_mode,'Color','red','LineStyle','--'); hold on;
semilogx(R,lambda,'Color','blue'); hold on;
%semilogx(R,lambda_theory_correct_N200,'Color','blue','LineStyle',':'); hold on;
xlabel('r')
ylabel('Lambda')
title('Asymmetry Modal Equation')
legend('GCM','QG','1D Equation')
savefig('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/figures/lambda_r_4pi.fig')
% figure(2); saveas(gcf,'Plots_Generals/Moist_Modal_Theory/Asymmetry_vs_r_Mode','epsc')

% figure(3)
% semilogx(r_factor,lambda_omega_mode,'Color','red','linewidth',1.6); hold on;
% %semilogx(r_factor,lambda_omegaQG_mode,'Color','red','LineStyle','--'); hold on;
% %semilogx(R,lambda,'Color','blue'); hold on;
% semilogx(R,lambda_theory_correct_N200,'Color','blue','linewidth',1.6); hold on;
% xlabel('Reduction factor r','interpreter','latex')
% ylabel('Asymmetry parameter $\lambda$','interpreter','latex')
% title('Asymmetry')
% legend('GCM','1-D Modal Theory')
% legend boxoff
% set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
% set(gca,'fontsize', 14);
% set(gca,'linewidth',1.5)
% ylim([0.5,1])
%saveas(gcf,'Plots_AOFD/Asymmetry_vs_r_ZGMode1','epsc')

% Make Video

close all;

v = VideoWriter('/net/aimsir/archive1/mkohl/Mode_Zurita_Gotor/videos/Mode_8pi.avi');
v.FrameRate = 1;
 
open(v);

for tt = 1:length(R)
    
    plot(x,w_field(:,tt)); 
    xlabel('x')
    ylabel('w')
    title(['r=',num2str(R(tt))]);
    xlim([0 L])
    
    
    frame = getframe(gcf);
    myVideo.FrameRate = 1;
    writeVideo(v,frame);
  
end
    

close(v);