function [asymmetry] = toy_model_2d_non_sine(r,kd)

% calculate the asymmetry based on the non-sinusoidal 2d toy model

kd2 = kd^2;
ku2 = (-1+sqrt(1+8*r.*kd2.*(1+2*kd2)))./(4*r); % alternative definition


asymmetry = 1./(1+kd2./ku2);



end

