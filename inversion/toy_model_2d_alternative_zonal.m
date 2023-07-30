function [w_final,asymmetry_parameter] = toy_model_2d_alternative_zonal(r,k)

% Calculate the Theoretical Prediction for the Asymmetry with dimensional 
% Wavenumber k as an input. 


% Define Parameters

%N = 1000; % number of gridpoints points
N = 300; % number of gridpoints points
Ny = N;

% Define Grid

p = 1; %length of periodic domain

% domain length L = 2*pi/k

x = linspace(0,2*p*pi/k,N); % define domain-length
y = linspace(0,2*p*pi/k,Ny); 

dx = 2*p*pi/(k*(N-1)); % define grid-spacing 
dy = 2*p*pi/(k*(Ny-1));

[X,Y] = meshgrid(x,y);


% Define Idealized RHS

%RHS = sin(k*X).*sin(k*Y); % RHS-forcing
RHS = sin(k*X);

% % periodicity - discard outer points
% 
RHS = RHS(1:end-1,1:end-1);

N = N-1;
Ny = Ny-1;

RHS = RHS(:); NN = length(RHS);

w_old = randn(length(RHS),1); % start with random guess

% Start the iterative-Inversion

E = 10; % Initialize Error-term
counter = 0;

while E>=1e-12
counter = counter + 1;
% Define the r-factor based on w
    
rr = ones(size(w_old)); rr(w_old>0) = r; 
 
% Define the A-matrix 
   
vp1x = sparse(rr(2:end)'); vp1x(N:N:end) = 0;
vm1x = sparse(rr(1:end-1)'); vm1x(N:N:end) = 0;

vl1x = sparse(zeros(1,NN-(N-1))); vl1x(1:N:end) = rr(N:N:end);
vr1x = sparse(zeros(1,NN-(N-1))); vr1x(1:N:end) = rr(1:N:end);

dA = sparse(diag(sparse(-2*rr'/dx^2-2*rr'/dy^2-1).*sparse(ones(1,NN))));
dAp1x = sparse(diag( sparse(1/dx^2*vp1x), 1 )); 
dAm1x = sparse(diag( sparse(1/dx^2.*vm1x), -1 ));
dAl1x = sparse(diag(1/dx^2*vl1x,(N-1)));
dAr1x = sparse(diag(1/dx^2*vr1x,-(N-1)));

dAp1y = sparse(diag(1/dy^2*sparse(rr(N+1:1:end)'), N ));
dAm1y = sparse(diag(1/dy^2*sparse(rr(1:1:NN-N)'),-N));


dAuy = sparse(diag(1/dy^2*sparse(rr(NN-(N-1):1:end)'),N*(Ny-1))); 
dAly = sparse(diag(1/dy^2*sparse(rr(1:1:N)'),-N*(Ny-1)));

A = (dA + dAp1x + dAm1x + dAl1x + dAr1x + dAp1y + dAm1y + dAuy + dAly);

% Carry out the inversion to find the new w

w_new = A\RHS;


%E = rms(w_new-w_old);
%E = max(w_new-w_old);

w_old = w_new;

rr = ones(size(w_old)); rr(w_old>0) = r; 
E = rms(Lap(N,Ny,dx,dy)*(rr.*w_old)-w_old-RHS);

% wtest = reshape(w_old,N,N); 
% plot(wtest(end/2,:)); hold on;

end
%counter
% toc

w = reshape(w_old,N,Ny);

% enforce periodicity

w_final = zeros(size(X));
w_final(1:end-1,1:end-1) = w;
w_final(end,:) = w_final(1,:);
w_final(:,end) = w_final(:,1);


asymmetry_parameter = Lambda(-w_final); % lambda form


end