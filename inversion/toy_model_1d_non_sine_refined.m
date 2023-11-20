function [asymmetry] = toy_model_1d_non_sine_refined(r,k)

% calculate the asymmetry based on the non-sinusoidal 2d toy model

% numerical solve

kd = k;

% % initial guess
% 
x0 = 1;
kuu = fsolve(@(x)root_1d_toy_model(x,kd,r),x0);
% res = kuu-kd*(1+kd^2)/(1+kuu^2*r)

% solve
%kuu = 1.3;
ku = fsolve(@(x) root_1d_toy_model_refined(x,r,k),kuu);
kd = k*ku/(2*ku-k);

res = ku-kd*(1+kd^2)/(1+ku^2*r)


asymmetry = 1/(1+kd/ku);

end