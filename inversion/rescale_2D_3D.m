function [rescale_sigma,rescale_b,rescale_sigma_av,rescale_b_av,...
 rescale_sigma_500,rescale_sigma_b_500] = rescale_2D_3D(u_eff,u_av,u_500,f_eff,delta_p,delta_p_500,S0,S0_mean,S0_500)

% Calculate the rescale factors 

day = 24*3600; km = 1000;

rescale_sigma = u_eff.*f_eff./(delta_p.*sqrt(-S0))*day;

rescale_b = sqrt(-S0).*delta_p./f_eff*1/km;

rescale_sigma_av = u_av.*f_eff./(delta_p.*sqrt(-S0_mean))*day;

rescale_b_av = sqrt(-S0_mean).*delta_p./f_eff*1/km;


rescale_sigma_500 = u_eff.*f_eff./(delta_p.*sqrt(S0_500/2))*day;

rescale_sigma_b_500 = sqrt(S0_500/2).*delta_p./f_eff*1/km;


%rescale_sigma_500 = u_500'.*f_eff./(delta_p_500*sqrt(S0_500/2))*day;

%rescale_sigma_b_500 = sqrt(S0_500/2).*delta_p_500./f_eff*1/km;


end

