function [F] = root_1d_toy_model_refined(x,r,k)

% solve the nonlinear constraint for x=ku

kd = k*x/(2*x-k);

F = x-kd*(1+kd^2)./(1+x^2*r);

end
