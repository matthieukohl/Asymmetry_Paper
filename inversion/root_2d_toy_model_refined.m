function [F] = root_2d_toy_model_refined(x,r,k)

% solve the nonlinear constraint for x=ku

kd = k*x/(2*x-k);

F = x^2 - kd^2*(1+2*kd^2)/(1+2*x^2*r);

end

