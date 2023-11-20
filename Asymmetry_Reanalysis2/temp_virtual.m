function [T_virtual] = temp_virtual(temp,p)
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

% density
T_virtual = temp.*(1.0+rs/gc_ratio)./(1.0+rs);



end