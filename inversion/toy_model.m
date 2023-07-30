function [w_final,asymmetry_parameter,asymmetry_parameter_2] = toy_model(R,k)

% Calculate the Theoretical Prediction for the Asymmetry with dimensional 
% Wavenumber k as an input. 


% Define Parameters

%N = 1000; % number of gridpoints points
N = 300; % number of gridpoints points

% Define Grid

p = 1; %length of periodic domain

% domain length L = 2*pi/k

x = linspace(0,2*p*pi/k,N); % define domain-length

dx = 2*p*pi/(k*(N-1)); % define grid-spacing 

% domain length L = 2*pi

%x = linspace(0,2*pi,N); % define domain-length

%dx = 2*pi/(N-1); % define grid-spacing 


% try inversions with fixed resolution and domain length
% adjust gridpoints accordingly

% p = 1;
% L = 2*p*pi/k;
% dx = 0.013;
% N = 1+L/dx;
% x = linspace(0,L,N);
% dx = x(2)-x(1);

% Define Idealized RHS

RHS = sin(k*x)'; % RHS-forcing

%Integral = trapz(x,RHS); % needed for BC

%RHS(1) = RHS(1) - 1/dx*Integral; % if Int RHS is non-zero

RHS = RHS(1:end-1); % periodicty; only need to solve for N-1 points. 

%RHS = RHS - mean(RHS);

% Pre-Initialize the r-factor, asymmetry parameter

r = zeros(length(RHS),1);

% Start the iterative-Inversion

E = 1; % Initialize Error-term

count = 0; % Initialize Counter
    
w_old = randn(length(RHS),1); % start with random guess
    
while E>=10e-12

count = count+1;    
    
% Define the r-factor based on w
    
 for j = 1:length(RHS)

    if w_old(j)>=0

        r(j) = R;
    else 
        r(j) = 1;
    end

 end
 
% Define the A-matrix 
    
dA = diag(-(2*r'/dx^2+1).*ones(1,length(RHS))); % diagonal matrix
dAp1 = diag( r(2:end)'/dx^2.*ones(1,length(RHS)-1), 1 ); % super-diagonal matrix
dAm1 = diag( r(1:end-1)'/dx^2.*ones(1,length(RHS)-1), -1 ); % sub-diagonal matrix
A = (dA + dAp1 + dAm1);
A(1,end) = r(end)/dx^2; %periodic BC
A(end,1) = r(1)/dx^2; %periodic BC

% Carry out the inversion to find the new w

w_new = A\RHS;

E = rms(w_new-w_old);

w_old = w_new;

end

w_final = [w_old;w_old(1)]; % enforce periodicity. 

asymmetry_parameter = asymmetry(w_final);

ind_up = w_final>0;
ind_down = w_final<0;

asymmetry_parameter_2 = 1-sum(ind_up(:))/(sum(ind_up(:))+sum(ind_down(:)));

% make spectral plot
% L = 2*pi/k;
% kk = 2*pi/L *[0:N/2, -N/2+1:-1];
% power = abs(fft(w_final)).^2;
% ind = N/2;
% 
% semilogx(kk(1:ind),power(1:ind));
% xlim([kk(1) kk(ind)]);

end