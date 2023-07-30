function [asymmetry] = toy_model_2d_non_sine_ku_input(r,ku)

% calculate the asymmetry based on the non-sinusoidal 2d toy model

ku2 = ku^2;


kd2 = 1/4*(-1 + sqrt(1+8*ku2+16*r*ku2^2));

asymmetry = 1./(1+kd2./ku2);



end