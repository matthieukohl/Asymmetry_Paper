function [r] = r_parameter(T,Level,level)
% Calculate the inferred r field using Fantini's Formula
% r = (theta)/(theta*)*(Gamma_M)/(Gamma_d)*(dtheta*/dp)/(dtheta/dp)

% Constants

Ra = 287.04; % Gas-constant Air
Rc = 461.92; % Gas-constant Water Vapour
cp = 1005.7; % Heat Capacity Air
cpc = 1850; % Heat Capacity Water Vapour
kappa = Ra/cp;
T_triple = 273.16; % Triple Point Temperature
p_triple = 611; % Triple Point Pressure
L = 2493*10^3; % Latent Heat Water Vapour
epsilon = 0.621; % Mass ratio Water Vapour / Air


% p_sat - saturation pressure water vapor

p_sat = p_triple*exp(-L/Rc*(1./T-1/T_triple));

% r_sat - saturated mixing ratio

r_sat = epsilon*p_sat./(Level-p_sat); 

% dry-lapse rate

Gamma_dry = Ra/cp;

% moist-lapse rate

Gamma_moist = Ra/cp*(1+L/(Ra*T).*r_sat)./(1+(cpc/cp+(L./(Rc*T)-1)*L./(cp*T)).*r_sat);

% potential temperature and gradient

theta = T.* (1e5 ./ Level) .^ kappa;

dTheta_dp = d_dp(theta,level);

% Saturation Equivalent Potential Temperature Bolton (1980) and gradient

theta_dl = T.*(1e5/(Level-p_sat)).^(kappa);

theta_star = theta_dl.*exp((3036./T-1.78).*r_sat.*(1+0.448*r_sat));

dTheta_star_dp = d_dp(theta_star,level);

% Calculation of r

r = (theta./theta_star).*(Gamma_moist./Gamma_dry).*(dTheta_star_dp./dTheta_dp);

r  = max(r(:,:,:),0);

r = min(r(:,:,:),1);

end

