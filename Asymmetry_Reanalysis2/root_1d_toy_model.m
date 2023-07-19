function [F] = root_1d_toy_model(x,kd,r)
%1d toy model: x = k_u

F = x-kd*(1+kd^2)./(1+x^2*r);

end

