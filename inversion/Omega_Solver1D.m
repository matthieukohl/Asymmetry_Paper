function [w,asymmetry] = Omega_Solver1D(RHS,r,N,dx)
% Solve 1-D moist omega equation (r(w)*w)_{xx} - w = RHS with equal grid

w_old = randn(length(RHS),1); % start with random guess

% Start the iterative-Inversion

E = 1; % Initialize Error-term
counter = 0;
while E>=10e-12
counter = counter + 1;

% Define the r-factor based on w
    
rr = ones(size(w_old)); rr(w_old>0) = r; 
 
% Define the A-matrix 
     
% dA = diag(-(2*rr'/dx^2+1).*ones(1,length(RHS))); % diagonal matrix
% dAp1 = diag( rr(2:end)'/dx^2.*ones(1,length(RHS)-1), 1 ); % super-diagonal matrix
% dAm1 = diag( rr(1:end-1)'/dx^2.*ones(1,length(RHS)-1), -1 ); % sub-diagonal matrix

dA = sparse(diag(sparse(-(2*rr'/dx^2+1)).*sparse(ones(1,length(RHS))))); % diagonal matrix
dAp1 = sparse(diag( sparse(rr(2:end)'/dx^2).*sparse(ones(1,length(RHS)-1)), 1 )); % super-diagonal matrix
dAm1 = sparse(diag( sparse(rr(1:end-1)'/dx^2).*sparse(ones(1,length(RHS)-1)), -1 )); % sub-diagonal matrix

A = (dA + dAp1 + dAm1);
A(1,end) = rr(end)/dx^2; %periodic BC
A(end,1) = rr(1)/dx^2; %periodic BC

% Carry out the inversion to find the new w

w_new = A\RHS;

E = rms(w_new-w_old);

w_old = w_new;

end

w = w_old;
asymmetry = Lambda(-w);
end

