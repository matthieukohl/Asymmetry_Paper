function [lambda] = toy_model_2d_non_sine_refined(r,k)

% solve the 2d toy model
kd = k;
kd2 = kd^2; 

% use the 2/k = 2/kd model as an initial guess

kuu2 = (-1+sqrt(1+8*r*kd2*(1+2*kd^2)))/(4*r);
kuu = sqrt(kuu2);

% solve

ku = fsolve(@(x) root_2d_toy_model_refined(x,r,k),kuu);
ku2 = ku^2;
kd = k*ku/(2*ku-k);
kd2 = kd^2;

% check if equation is solved
res = ku2-kd2*(1+2*kd2)/(1+2*ku2*r)

lambda = 1/(1+kd2/ku2); 

end

