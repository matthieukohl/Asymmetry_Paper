function [asymmetry] = toy_model_1d_non_sine(r,kd)

% calculate the asymmetry based on the non-sinusoidal 2d toy model

% numerical solve

x0 = 1;
ku = fsolve(@(x)root_1d_toy_model(x,kd,r),x0);

% analytic version based on wolfram alpha results

%express = nthroot(27*r^2*kd^3+27*r^2*kd+sqrt(108*r^3+(27*r^2*kd^3+27*r^2*kd)^2),3);

%express = (27*r^2*kd^3+27*r^2*kd+sqrt(108*r^3+(27*r^2*kd^3+27*r^2*kd)^2))^(1/3);

%ku = express/(3*r*2^(1/3)) - 2^(1/3)/express;

% ku2 = (1+1i*sqrt(3))/(2^(2/3)*express)...
%      -(1-1i*sqrt(3))*express/(6*2^(1/3)*r);
% 
% ku3 = (1-1i*sqrt(3))/(2^(2/3)*express)...
%      -(1+1i*sqrt(3))*express/(6*2^(1/3)*r);
% 
% 
% analytic version based on complete cubic root formula 
% (Wolfram: https://mathworld.wolfram.com/CubicFormula.html)
% also useful: https://www.e-education.psu.edu/png520/m11_p6.html

% a2 = 0; a1 = 1/r; a0 = -(kd+kd^3)/r;
% %a2 = 3/2; a1 = -11/2; a0 = -3; test purely real solutions
% 
% Q = (3*a1-a2^2)/9;
% R = (9*a2*a1-27*a0-2*a2^3)/54;
% 
% D = Q^3 + R^2;
% 
% % S = (R+sqrt(D))^(1/3);
% % T = (R-sqrt(D))^(1/3);
% 
% % beware the cubic roots of negative numbers;
% % nthroot picks the real one
% 
% S = nthroot(R+sqrt(D),3); 
% T = nthroot(R-sqrt(D),3); 
% 
% ku = -1/3*a2 + S + T;
% % kuu2 = -1/3*a2 - 1/2*(S+T) + 1/2*1i*sqrt(3)*(S-T);
% % kuu3 = -1/3*a2-1/2*(S+T)-1/2*1i*sqrt(3)*(S-T);
% 
% 
asymmetry = 1/(1+kd/ku);

end