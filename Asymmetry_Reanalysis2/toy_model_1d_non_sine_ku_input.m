function [asymmetry] = toy_model_1d_non_sine_ku_input(r,ku)

% calculate the asymmetry based on the non-sinusoidal 2d toy model

% % numerical solve
% 
x0 = 1;
kd = fsolve(@(x)root_1d_toy_model_ku_input(x,ku,r),x0);

asymmetry = 1/(1+kd/ku);

end