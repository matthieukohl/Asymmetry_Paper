function F = root_finder(k,r,x)

% root find the transition point b for the toy-model
% imposing continuity of (rw)_x as boundary condition

cd = -sin(k*x(1))/(1+k^2)*1/cosh(x(1)+pi/(2*k));
cu = -sin(k*x(1))/(1+r*k^2)*1/cosh(1/sqrt(r)*(x(1)-pi/(2*k)));

xb = pi/k+x(1);


wdx = -k*cos(k*xb)/(1+k^2)+cd*sinh(xb-pi/(2*k));

wux = -k*cos(k*xb)./(1+r*k^2)+cu*1/sqrt(r)*sinh(1/sqrt(r)*(xb-3*pi/(2*k)));



F = r*wux-wdx;


end

