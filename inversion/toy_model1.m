function [asymmetry_parameter] = toy_model1(r,k)

% define physical constants

f = 1.0491e-04; % s^(-1) Coriolis Parameter

N2 = 1.1870e-04; % s^(-2) Static Stability

H = 10; % km Vertical Scale

L=sqrt(N2)*H/f;  %Deformation Horizontal Scale

%define model parameters

a = 4; % height scaling parameter

k = k*L/(6371*cosd(46)); % forcing Wavenumber

% Define Parameters

R = r; % reduction factor r

N = 300; % number of gridpoints points

% Define Grid

p = 1; %length of periodic domain

x = linspace(0,2*p*pi/k,N); % define domain-length

dx = 2*p*pi/(k*(N-1)); % define grid-spacing 

% Define RHS

RHS = sin(k*x)'; % RHS-forcing

Integral = trapz(x,RHS); % needed for BC

RHS(1) = RHS(1) - 1/dx*Integral; % if Int RHS is non-zero

RHS = RHS(1:end-1); % periodicty; only need to solve for N-1 points. 

% Pre-Initialize the r-factor, asymmetry parameter

r = zeros(length(RHS),1);

% Start the iterative-Inversion

E = 1; % Initialize Error-term

count = 0; % Initialize Counter
    
w_old = randn(length(RHS),1); % start with random guess
    
while E>=10e-6

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
    
dA = diag(-(2*r'/dx^2+a).*ones(1,length(RHS))); % diagonal matrix
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

end