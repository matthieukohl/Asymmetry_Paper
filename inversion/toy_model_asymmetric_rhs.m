function [w_final,asymmetry_parameter] = toy_model_asymmetric_rhs(R,k)

% Calculate the Theoretical Prediction for the Asymmetry with dimensional 
% Wavenumber k as an input. 

% Define Parameters

N = 300; % number of gridpoints points

% Define Grid

p = 1; %length of periodic domain

x = linspace(0,2*p*pi/k,N); % define domain-length

dx = 2*p*pi/(k*(N-1)); % define grid-spacing 

% Define Idealized RHS with slight asymmetry

RHS = sin(k*x)'; % RHS-forcing

%RHS(RHS>0) = 0.01*RHS(RHS>0);
RHS(RHS<0) = 5*RHS(RHS<0);

%RHS = pearsrnd(0,1,-1,3,length(x),1);

%Integral = trapz(x,RHS); % needed for BC

%RHS(1) = RHS(1) - 1/dx*Integral; % if Int RHS is non-zero

RHS = RHS(1:end-1); % periodicty; only need to solve for N-1 points. 

RHS = RHS - mean(RHS);

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

% remove the mean because now RHS no longer has mean zero

w_final = w_final - mean(w_final);

asymmetry_parameter = asymmetry(w_final);


end