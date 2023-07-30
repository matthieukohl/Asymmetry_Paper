function [omega_adv, omega_diab] = inversion_toy_rhs(u,v,omega,T,r,sigma,f,beta,level,Phi,phi,dphi,lambda,dlambda)
% Invert the moist-omega equation to find w due to advection vs. 
% w due to advection 

% Define Constants

R = 6371000.0; % Planetary Radius 
Ra = 287.04; % Gas-constant Air

% Define domain sizes

Nlat = size(omega,1);
Nlon = size(omega,2);
Nlevel = size(level,1);

% Predefine Terms

[A, B, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi] = deal(zeros(Nlat,Nlon,Nlevel));
    
[rhs, omega_b] = deal(zeros(Nlat,Nlon,Nlevel));
[omega_adv,omega_diab] = deal(nan(Nlat,Nlon,Nlevel));


% term A (Q vector)

for k = 1 : length(level)
    
    du_dlambda(:, :, k) = d_dlambda(u(:, :, k), dlambda);
    dv_dlambda(:, :, k) = d_dlambda(v(:, :, k), dlambda);
    du_dphi(:, :, k) = d_dphi(u(:, :, k), dphi);
    dv_dphi(:, :, k) = d_dphi(v(:, :, k), dphi);

    dT_dlambda(:, :, k) = d_dlambda(T(:, :, k), dlambda);
    dT_dphi(:, :, k) = d_dphi(T(:, :, k), dphi);
    Qx(:, :, k) = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda(:, :, k) .* dT_dlambda(:, :, k) + ...
        (1 ./ cos(Phi))    .* dv_dlambda(:, :, k) .* dT_dphi(:, :, k));
    Qy(:, :, k) = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi(:, :, k)    .* dT_dlambda(:, :, k) + ...
                              dv_dphi(:, :, k)    .* dT_dphi(:, :, k));
    div_temp = div(Qx(:, :, k), Qy(:, :, k), Phi, dphi, dlambda);
    A(:, :, k) = 2 * Ra / level(k) * div_temp;
end

clear('temp_x', 'temp_y', 'div_temp');

% term B (planetary vorticity)

temp = d_dp(v(:, :, :), level);
for j = 1 : Nlat
    B(j, :, :) = f(j) * beta(j) * temp(j, :, :);
end
clear('temp');

% right hand side

rhs(:, :, :) = A(:, :, :) + B(:, :, :); 

% modify rhs to have sinusoidal rhs profile 
% rhs = rhs500(x,y)*sin(pi*p/p_surf)

% rhs_500 = rhs(:,:,11);
% 
% RHS_500 = repmat(rhs_500,1,1,length(level));
% 
% Level = repmat(reshape(level, 1, 1, length(level)), length(Nlat), length(Nlon), 1);
% 
% rhs = RHS_500.*sin(pi*Level/level(1));

% modify rhs to have constant profile

rhs_500 = rhs(:,:,11);

RHS_500 = repmat(rhs_500,1,1,length(level));

rhs = RHS_500;



% cut to lowest non-NaN level in sigma

for bound=1:length(level)

b = isnan(sigma(bound)); 

if any(b(:))==0
    break
end 

end

rhs_place = rhs(:,:,bound:end);
omega_b_place = omega_b(:,:,bound:end);
sigma_place = sigma(bound:end);
level_place = level(bound:end);
r_place = r(:,:,bound:end);
omega_place = omega(:,:,bound:end);

% Start Iterative Inversion with Random field

count = 0;

w_old = randn(Nlat,Nlon,length(level_place)); 

w_new = SIP_diagnostic1(R, f, Nlon, Nlat, length(level_place), phi, lambda, level_place, ...
sigma_place(:), dphi, dlambda, rhs_place(:,:,:),  omega_b_place(:,:,:));

omega_adv(:,:,bound:end) = w_new;

%w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%sigma_place(:), dphi, dlambda, rhs_place(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));

diff = w_new - w_old;

% Pre-Allocation

sigma_r = deal(zeros(Nlat,Nlon,length(level_place)));
del2_J = deal(zeros(Nlat,Nlon,length(level_place)));

while rms(diff(:))>=10^(-4) && count<350
    
count = count + 1;
       
 for l=1:Nlat
    for m=1:Nlon
      for n= 1:length(level_place)
          if w_new(l,m,n)>=0
           sigma_r(l, m, n) = 0*sigma_place(n); 
          else 
           sigma_r(l, m, n) = (1-r_place(l,m,n))*sigma_place(n);  
          end 
      end
    end
 end

J = w_new.*sigma_r;


for k=1:size(level_place)
    
    del2_J(:,:,k) = spherical_laplacian(J(:,:,k),phi,dphi,dlambda);
    
end
 
rhs_mod = rhs_place+del2_J; 

w_old = w_new;

w_new = SIP_diagnostic1(R, f, Nlon, Nlat, length(level_place), phi, lambda, level_place, ...
                   sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_b_place(:,:,:));

%w_new = SIP_diagnostic2(R, f, length(lon), length(phi), length(level_place), phi, lambda, level_place, ...
%sigma_place(:), dphi, dlambda, rhs_mod(:,:,:), omega_place(:,:,:),omega_b_place(:,:,:));                            
                             
                             
diff = w_new - w_old; 
rms(diff(:));
error(count) = rms(diff(:));
asym(count) = Lambda(w_new(:,:,1:8));
end
count
omega_diab(:,:,bound:end) = w_new;

end
