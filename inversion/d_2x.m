function [A] = d_2x(N,dx)

% Calculate the second derivative of f

a = -2/dx^2; a1 = 1/dx^2; 

dA = diag((a).*ones(1,N)); % diagonal entries
dAp1 = diag((a1).*ones(1,N-1),1); % i+1 entries
dAm1 = diag((a1).*ones(1,N-1),-1); % i-1 entries

A = (dA + dAp1 + dAm1); 

A(1,end) = a1; %periodic BC
A(end,1) = a1; %periodic BC


end

