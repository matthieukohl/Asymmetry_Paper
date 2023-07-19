function [r] = r_parameter1(temp,p,level)
% Calculate r=1-dTheta/dp_theta* / dTheta/dp based on O'Gorman formulation. 

% constants
Rd       = 287.04;      % gas constant for dry air [J/kg/K]
Rv       = 461.5;       % gas constant for water vapor [J/kg/K]
cpd      = 1005.7;      % specific heat dry air [J/kg/K]
cpv      = 1870;        % specific heat water vapor [J/kg/K]
g        = 9.80665;     % gravitational acceleration [m/s^2]
p0       = 1e5;         % reference pressure [Pa]
kappa    = Rd/cpd;
gc_ratio = Rd/Rv;

% saturation vapor pressure [Pa]
Tc = temp-273.15;
es = 611.20*exp(17.67*Tc./(Tc+243.5)); % Bolton 1980 equation 10

% latent heat of condensation [J/kg]
L = (2.501-0.00237*Tc)*1e6; % Bolton 1980 equation 2

% saturation mixing ratio
rs = gc_ratio*es./(p-es);

% saturation specific humidity 
qs = rs./(1+rs);

% potential temperature
exponent = kappa.*(1+rs./gc_ratio)./(1+rs.*cpv./cpd);
theta    = temp.*(p0./p).^exponent;

% derivative of potential temperature with respect to pressure
dtheta_dp = d_dp(theta, level);

% density
temp_virtual = temp.*(1.0+rs/gc_ratio)./(1.0+rs);
rho          = p/Rd./temp_virtual;

% moist adiabatic lapse rate
malr = g/cpd*(1+rs)./(1+cpv/cpd.*rs) ...
    .*(1+L.*rs./Rd./temp)...
    ./(1+L.^2.*rs.*(1+rs/gc_ratio)./(Rv.*temp.^2.*(cpd+rs*cpv)));

% derivative of potential temperature wrt pressure along a moist adiabat
% (neglects small contribution from vertical variations of exponent)
dtemp_dp_ma  = malr/g./rho;
dtheta_dp_ma = dtemp_dp_ma.*theta./temp-exponent.*theta./p;

% Calculation of r-parameter

r = 1-dtheta_dp_ma./dtheta_dp;

% set negative values to zero
r = max(r,0);

%r(r<0) = 1;

% set values larger than 1 to 1
r = min(r,1);


end

