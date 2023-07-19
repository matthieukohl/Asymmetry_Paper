function [w_final,asymmetry_parameter] = toy_model_nondimensional_generated_spectrum(R,RHS,dx)

% Calculate the Theoretical Prediction for the Asymmetry with dimensional 
% Wavenumber k as an input. 

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
    
% N = length(RHS); 
% r = ones(N,1);
% r(w_new>0) = R;
% E = d_2x(N,dx)*(r.*w_new)-w_new-RHS;
% E = E/max(abs(w_new));
% E = rms(E);

w_old = w_new;

end

w_final = [w_old;w_old(1)]; % enforce periodicity. 

asymmetry_parameter = asymmetry(w_final);

end