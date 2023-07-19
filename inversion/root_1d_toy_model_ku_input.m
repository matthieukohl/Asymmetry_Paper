function [F] = root_1d_toy_model_ku_input(x,ku,r)
%1d toy model: x = k_d

F = x - ku*(1+r*ku^2)/(1+x^2);

end