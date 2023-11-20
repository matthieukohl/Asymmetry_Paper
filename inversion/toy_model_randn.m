function [asymmetry_parameter] = toy_model_randn(sigma,r,lat,lon,level,phi,dphi,lambda,dlambda,strat_level)

% Calculate lambda using a randn RHS.

% Define Constants
R = 6371000.0;
Omega = 7.2921e-5; % Rotation Rate
f0 = 2 * Omega * sind(lat(round(length(lat)/2,0))); % Coriolis Parameter
runs = 20;

% Initialize lambda_run

lambda_run = zeros(length(runs),1);

for tt = 1:runs
    
% Initialize Gauss-Random RHS
[sigma_r,omega_b,del2_J] = deal(zeros(length(lat),length(lon),length(level)));

rhs = randn(length(lat),length(lon),length(level));

% w = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
% sigma(:), dphi, dlambda, rhs(:,:,:), omega_b(:,:,:));

w = randn(length(lat),length(lon),length(level));

% Start Iterative Inversion

for s=1:30
       
 for l=1:length(lat)
    for m=1:length(lon)
      for n= 1:length(level)
          if w(l,m,n)>=0
           sigma_r(l, m, n) = 0*sigma(n); 
          else 
           sigma_r(l, m, n) = (1-r(l,m,n))*sigma(n);  
          end 
      end
    end
 end

J = w.*sigma_r;


for k=1:size(level)
    
del2_J(:,:,k) = spherical_laplacian(J(:,:,k),phi,dphi,dlambda);
    
end
 
rhs_mod = rhs+del2_J;


w = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
             sigma(:), dphi, dlambda, rhs_mod(:,:,:), omega_b(:,:,:));


end


lambda_run(tt) = Lambda(w(:,:,1:strat_level,:));

end

asymmetry_parameter = mean(lambda_run(tt));

end

