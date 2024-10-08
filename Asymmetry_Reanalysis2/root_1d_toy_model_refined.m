function [F] = root_1d_toy_model_refined(x,r,k)

%1d toy model: x = k_u

% solve the nonlinear constraint for x=ku

kd = k*x/(2*x-k);

F = x-kd*(1+kd^2)/(1+x^2*r);

%F = x*(1+r*x^2) - kd*(1+kd^2);

% test both sides

% k = 1;
% x = linspace(0.7,5,100);
% kd = k*x./(2*x-k);
% r = 0.01;
% test = kd.*(1+kd.^2)./(1+x.^2*r);

%plot(x,x); hold on; plot(x,test);

end