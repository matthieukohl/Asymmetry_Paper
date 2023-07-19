function [lambda,w_final] = mode_asymmetry(R,N,L)

%% Tilted 2-Layer QG-Model with Latent Heating and no Internal PV-Gradients
% Tilted Interfaces (to produce zero internal PV-gradients) h_b = h_t = h = y
% and linearization around a mean shear state tau = -y produces the modelling
% equation: 
% (i) d_t phi_xx + tau_xxx - h*tau_x = 0
% (ii) d_t tau_xx + phi_xxx - h*phi_x + w = 0
% (iii) (r(w)w)_xx - w = 2*phi_xxx - h*phi_x
% h = 0 for the classic untilted model; h = 1 for the tilted model

%tic
% add path

addpath('/net/halo/disk28/disk7/mkohl/Mode_Kohl_OGorman');

% Define Grid

%N = 400;
%L = 16*pi;
x = linspace(0,L,N+1);
dx = L/(N);


% Define time-step

t_final = 200;
% tN = 16000;
% dt = t_final/tN;


% untilted h = 0, tilted h = 1;

h = 0;

% Initial Conditions

IC = 'rand';
%IC = 'sinusoidal';

% define variables

sigma = zeros(length(R),1);  

w_final = zeros(N+1,length(R));


% Rescale the Fields

options = odeset('Events',@myEvent);


for ii = length(R):length(R)



% Initialize Values

% Initial Conditions: Choose between random and sinusoidal guess

if strcmp(IC,'rand') == 1
% random-guess
phi =   randn(1,N); phi = phi - mean(phi); 
tau =   randn(1,N); tau = tau - mean(tau); 
else

% sinusoidal-guess

%phi = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); phi = phi(1:end-1);
%tau = sin(x)+sin(2*x)+sin(4*x)+sin(6*x)+sin(8*x)+sin(10*x); tau = tau(1:end-1);

phi = sin(x); phi = phi(1:end-1);
tau = sin(x); tau = tau(1:end-1);

end

Phi0 = [tau,phi]';

% Define inverse matrix

A2 = d_2x(N,dx); A2(1,:) = 1; A2inv = inv(A2);

%Propagator: 
total_time = 0;
g = [];
W = [];
s_mode = [];
skew_rhs = [];
time = 0;
nrescale = 0;

while total_time < t_final
[t,y] = ode45(@(t,y) Prop(t,y,A2inv,R(ii),h,N,dx),[0 t_final-total_time],Phi0,options);
Phi0 = y(end,:)/100;
total_time = t(end)+total_time;
g = cat(1,g,y);

% calculate the growthrate over the time of a rescaling event
% s = 1/(t(end)-t(1))*log(rms(y(end,:))/rms(y(1,:)));
% 
% s_mode = cat(1,s_mode,s);

% calculate the growthrate for every timestep

time = cat(1,time,t(2:end)+time(end));

s = zeros(length(t),1);
for kk = 2:length(t)
    s(kk) = 1/(t(kk)-t(kk-1))*log(rms(y(kk,:))/rms(y(kk-1,:)));

end
s_mode = cat(1,s_mode,s(2:end));

% % calculate the skewness of the RHS for h = 0
% 
% for jj = 2:size(y,1)
% RHS = 2*d_1x(N,dx)*y(jj,N+1:end)';
% skew_rhs = cat(1,skew_rhs,skewness(-RHS));
% end

nrescale = nrescale +1;
end

dt = t(end)-t(end-1);

Phi = y(:,N+1:end)'; Tau = y(:,1:N)';

% Last step w-inversion
Phi_end = Phi(:,end); Phi_end(1) = 0;
Phi_p_end = Phi(:,end-1); Phi_p_end(1) = 0;

Tau_end = Tau(:,end); Tau_end(1) = 0;
Tau_p_end = Tau(:,end-1); Tau_p_end(1) = 0;

phi_end = A2inv*Phi_end; phi_p_end = A2inv*Phi_p_end;
tau_end = A2inv*Tau_end; tau_p_end = A2inv*Tau_p_end;

RHS = 2*d_1x(N,dx)*Phi(:,end)-h*d_1x(N,dx)*phi_end;
RHS1 = 2*d_1x(N,dx)*Phi(:,end);
RHS2 = -d_1x(N,dx)*phi_end;

w = Omega_Solver(RHS,R(ii),N,dx);
r = r_factor(w,R(ii));

% final fields 

w_final(:,ii) = [w;w(1)];

% growthrate: average the growth rate over the last five rescale events
% sigma(ii) = mean(s_mode(end-4:end));

% growthrate: average the stepwise growthrate over the last 5 simulation
% steps
[val,index] = min(abs(time-195));
sigma(ii) = mean(s_mode(index:end));

% Calculate Asymmetry

lambda = asymmetry(w_final(:,ii));


end
%toc



end

