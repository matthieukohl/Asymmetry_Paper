function [F] = root_2d_toy_model_refined(x,r,k)

% solve the nonlinear constraint for x=ku

kd = k*x/(2*x-k);

F = x^2 - kd^2*(1+2*kd^2)/(1+2*x^2*r);


% % test cross over
% k = 1;
% x = linspace(0.6,5,100);
% kd = k*x./(2*x-k);
% r = 0.01;
% test = kd.*(1+kd.^2)./(1+x.^2*r);
% test2 = kd.^2.*(1+2*kd.^2)./(1+2*x.^2*r);
% 
% %plot(x,x); hold on; plot(x,test);
% 
% plot(x,x.^2); hold on; plot(x,test2);


end

