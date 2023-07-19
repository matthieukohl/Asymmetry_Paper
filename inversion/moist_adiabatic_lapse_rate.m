function lapse_rate = moist_adiabatic_lapse_rate(T, p, type)
% MOIST_ADIABATIC_LAPSE_RATE Returns saturated moist-adiabatic lapse rate.
%  
% Units are K / m.

  g               = PARS('gravity');
  cpd             = PARS('cp');
  cpv             = PARS('cp_v');
  Rd              = PARS('gas_constant');
  Rv              = PARS('gas_constant_v');
  gc_ratio        = Rd / Rv;
  
  persistent not_first_call
  
  [es, qs, rs, L] = saturation_thermodynamics(T, p, type);

    
  switch lower(type)
   case 'simple'
    if isempty(not_first_call)
      disp(['MOIST_ADIABATIC_LAPSE_RATE: Using simplified moist adiabatic ' ...
            'lapse rate [as in dargan_betts_miller]'])
      not_first_call = 1;
    end

    % cf. Holton (p. 503) [used in dargan_betts_miller], obtained
    % by using rs<<1 in the more accurate expression below (though
    % this approximation is not necessary)
    lapse_rate      = g/cpd ...
        .* (1 + L.*rs ./ Rd ./ T)...
        ./ (1 + L.^2 .* rs ./(cpd * Rv .* T.^2));
    
   otherwise
    % cf. Emanuel (p. 131)
    lapse_rate      = g/cpd * (1 + rs) ./ (1 + cpv/cpd.*rs) ...
        .* (1 + L.*rs ./ Rd ./ T)...
        ./ (1 + L.^2 .* rs .* (1 + rs/gc_ratio)./(Rv .* T.^2 .* (cpd + rs*cpv)));
  end
     
  
 
