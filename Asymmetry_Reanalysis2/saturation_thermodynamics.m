function [es, qs, rs, L] = saturation_thermodynamics(T, p, type)
% SATURATION_THERMODYNAMICS  Saturation quantities for water.
%
% [es, qs, rs, L] = SATURATION_THERMODYNAMICS(TEMP, PRESS, type) returns
% the effective saturation vapor pressure, specific humidity, mixing
% ratio, and latent heat of vaporization/sublimation. Depending on
% type, is uses different formulations:
%    'era': Formulation used by ECMWF based on the modified
%           Tetens formula.
%  
% 'simple': Formulation used in idealized GCM experiments by Paul
%           O'Gorman and others with fixed latent heat. 
  
  Rd          = PARS('gas_constant');
  Rv          = PARS('gas_constant_v');
  gc_ratio    = Rd/Rv;   % ratio of gas constants for dry air and water vapor

  if nargin < 3
    type = 'simple';
  end

  persistent been_verbose
  if isempty(been_verbose)
    disp(['SATURATION_THERMODYNAMICS: Using ', upper(type), ...
          ' formulation for computation of saturation quantities.'])
    been_verbose = 1;
  end
    
  switch lower(type)
    
   case 'era'

    % ECMWF formulation described by Simmons et al. (1999: QJRMS, 125,
    % 353--386), which uses saturation over ice for temperatures
    % less than 250 K and a quadratic interpolation between
    % saturation over ice and over liquid water for temperatures
    % between 250 K and 273 K
    
    % coefficients in Tetens saturation vapor pressure es = es0 * exp(a3 * (T-T0)/(T-a4))      
    es0       = 611.21;  % saturation vapor pressure at T0 (Pa)
    T0        = 273.16;  % (K)
    Ti        = T0 - 23; % (K)

    a3l       = 17.502;  % liquid water (Buck 1981)
    a4l       = 32.19;   % (K)
  
    a3i       = 22.587;  % ice (Alduchov and Eskridge 1996)
    a4i       = -0.7;    % (K)
  
    % saturation vapor pressure over liquid and ice
    esl       = es0 .* exp(a3l .* (T - T0)./(T - a4l));
    esi       = es0 .* exp(a3i .* (T - T0)./(T - a4i));

  
    % latent heat of sublimation 
    Ls0       = 2.834e6;  % latent heat of sublimation  (J / kg) [+- 0.01 error for 173 K < T < 273 K]
    Ls        = Ls0 * ones(size(T));

    % latent heat of vaporization
    Lv0       = 2.501e6;  % latent heat of vaporization at triple point (J / kg)
    cpl       = 4190;     % heat capacity of liquid water (J / kg / K)
    cpv       = 1870;     % heat capacity of water vapor (J / kg / K)			  
    Lv        = Lv0 - (cpl - cpv) .* (T - T0);
  
    % compute saturation vapor pressure and latent heat over liquid/ice mixture
    iice      = find(T <= Ti);
    iliquid   = find(T >= T0);
    imixed    = find(T > Ti & T < T0);
  
    es        = ones(size(T))*NaN;
    L         = ones(size(T))*NaN;
    if ~isempty(iice)
      es(iice)  = esi(iice);
      L(iice)   = Ls(iice);
    end
  
    if ~isempty(iliquid)
      es(iliquid) = esl(iliquid);
      L(iliquid)  = Lv(iliquid);
    end
    
    if ~isempty(imixed)
      a         = ( (T(imixed) - Ti)./(T0 - Ti) ).^2;
      es(imixed)= (1-a) .* esi(imixed) + a .* esl(imixed);
      L(imixed) = (1-a) .* Ls(imixed) + a .* Lv(imixed);
    end

   case 'simple'
    
    % assuming constant latent heat of vaporization and only one
    % vapor-liquid phase transition
    T0        = 273.16;
    es0       = 610.78;
    L         = PARS('latent_heat_v');   
    
    % saturation vapor pressure
    es        = es0 * exp(L/Rv*(1.0/T0 - 1.0./T));
    
   otherwise
    error(['Unknown type, ', type, 'of computing saturation quantities.'])
  end
    
  % saturation mixing ratio
  rs          = gc_ratio * es ./ (p - es);
  
  % saturation specific humidity 
  qs          = rs ./ (1 + rs); 
  
  
